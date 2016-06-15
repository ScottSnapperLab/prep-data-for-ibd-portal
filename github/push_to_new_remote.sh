#!/bin/bash



github-create () {

    repo_name=$1

    dir_name=$(basename $(pwd))

    CURL=$(which curl)

    [ -d ".git" ] && INIT_GIT="y" || INIT_GIT="n"



    if [[ $CURL == *anaconda* ]]; then
        echo "[git-create][error]: You appear to be using an Anaconda-based version of curl; exit all anaconda environments and start again."
        invalid_credentials=1
    fi



    if [ "$repo_name" = "" ]; then
        echo "[git-create][info]: Repo name (hit enter to use '$dir_name')?"
        read repo_name
    fi

    if [ "$repo_name" = "" ]; then
        repo_name=$dir_name
    fi

    username=$(git config github.user)
    if [ "$username" = "" ]; then
        echo "[git-create][error]: Could not find username; run 'git config --global github.user <username>'"
        invalid_credentials=1
    fi

    token=$(git config github.token)
    if [ "$token" = "" ]; then
        echo "[git-create][error]: Could not find token; run 'git config --global github.token <token>'"
        invalid_credentials=1
    fi

    if [[ $invalid_credentials == "1" ]]; then
        return 1
    fi

    if [[ $INIT_GIT == "n" ]]; then
        echo "[git-create][info]: No git init detected; creating and performing intial commit..."
        git init
        git add .
        git commit -m "First commit"
    fi



    echo "[git-create][info]: Creating Github repository '$repo_name' ..."
    curl -u "$username:$token" https://api.github.com/user/repos -d '{"name":"'$repo_name'"}' # > /dev/null 2>&1
    echo "[git-create][info]:  done."

    echo "[git-create][info]: Pushing local code to remote ..."
    git remote add origin git@github.com:$username/$repo_name.git # > /dev/null 2>&1
    git push -u origin master # > /dev/null 2>&1
    echo "[git-create][info]:  done."
}




github-create "prep-data-for-ibd-portal"
