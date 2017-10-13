#!/bin/bash

set -e
set -x

if [[ "$(uname -s)" == 'Darwin' ]]; then

    # Install some custom requirements on OS X
    # e.g. brew install pyenv-virtualenv
    sw_vers
    brew update || brew update
    brew outdated openssl || brew upgrade openssl
    brew install openssl@1.1
    git clone --depth 1 https://github.com/pyenv/pyenv ~/.pyenv
    PYENV_ROOT="$HOME/.pyenv"
    PATH="$PYENV_ROOT/bin:$PATH"
    eval "$(pyenv init -)"
    case "${TOXENV}" in
        py34)
            # Install some custom Python 3.4 requirements on OS X
            pyenv install 3.4.6
            pyenv global 3.4.6
            ;;
        py35)
            # Install some custom Python 3.5 requirements on OS X
            pyenv install 3.5.3
            pyenv global 3.5.3
            ;;
        py36)
            # Install some custom Python 3.6 requirements on OS X
            pyenv install 3.6.1
            pyenv global 3.6.1
            ;;
    esac
    pyenv rehash
fi
