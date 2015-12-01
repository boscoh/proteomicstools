DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $DIR

BASENAME=${0##*/}
BASE=${BASENAME%.*}

python ${BASE}.py -i

if [ "$(uname)" == "Darwin" ]; then
    echo -n -e "]0;${BASE}-TKFORM"
    osascript -e 'tell application "Terminal" to close (every window whose name contains "'${BASE}'-TKFORM")' &
fi
