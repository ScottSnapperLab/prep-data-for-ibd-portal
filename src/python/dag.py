#!/usr/bin/env python
"""Provide code for manipulating DAG and rule graphs."""

# Imports
import textwrap
from collections import OrderedDict

# Metadata
__author__ = "Gus Dunn"
__email__ = "w.gus.dunn@gmail.com"


# DAG and rulegraph stuff
def digest_node_line(line):
    """Return OrderedDict of relevant line parts."""
    l = line.strip()

    d = OrderedDict()
    d["num"], fields = l.split('[')
    fields = fields.replace('rounded,dashed','rounded-dashed')
    fields = fields.rstrip('];').split(',')
    fields[-1] = fields[-1].replace('rounded-dashed','rounded,dashed')
    for f in fields:
        key, value = f.split('=')
        d[key.strip()] = value.strip().replace('"','').replace("'","")

    return d

def should_ignore_line(line, strings_to_ignore):
    """Return true if line contains a rule name in `rule_names`."""
    for string in strings_to_ignore:
        if string in line:
            return True

    return False

def recode_graph(dot, new_dot, pretty_names, rules_to_drop, color=None, use_pretty_names=True):
    """Change `dot` label info to pretty_names and alter styling."""
    if color is None:
        color = "#50D0FF"

    node_patterns_to_drop = []

    with open(dot, mode='r') as dot:
        with open(new_dot, mode='w') as new_dot:
            for line in dot:
                if '[label = "' in line:

                    # Add pretty names and single color IF pretty names are provided.
                    data = digest_node_line(line=line)
                    rule_name = data['label']

                    if use_pretty_names:
                        pretty_name = textwrap.fill(pretty_names[rule_name], width=40).replace('\n', '\\n')
                        full_name = "[{rule_name}]\\n{pretty_name}".format(rule_name=rule_name,pretty_name=pretty_name)
                        data['label'] = full_name
                        data['color'] = color
                    else:
                        pass

                    fields = ', '.join(['{k} = "{v}"'.format(k=k,v=v) for k, v in data.items()][1:])

                    if should_ignore_line(line, strings_to_ignore=rules_to_drop):
                        node_patterns_to_drop.append("\t{num} ->".format(num=data['num']))
                        node_patterns_to_drop.append("-> {num}\n".format(num=data['num']))
                        continue

                    new_line = """\t{num}[{fields}];\n""".format(num=data['num'],fields=fields)

                    new_dot.write(new_line)
                else:
                    if should_ignore_line(line, strings_to_ignore=node_patterns_to_drop):
                        continue
                    elif "fontname=sans" in line:
                        line = line.replace("fontname=sans","fontname=Cantarell")
                        line = line.replace("fontsize=10","fontsize=11")
                        new_dot.write(line)
                    else:
                        new_dot.write(line)
