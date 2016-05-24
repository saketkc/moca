#!/bin/bash
set -euo pipefail
cd ~
wget http://meme-suite.org/meme-software/4.10.2/meme_4.10.2.tar.gz
tar -zxvf meme_4.10.2.tar.gz 1>/dev/null 2>/dev/null
cd meme_4.10.2
./configure --prefix="$HOME/meme_bin" 1>/dev/null 2>/dev/null
make install 1>/dev/null 2>/dev/null
cd ~
wget http://meme-suite.org/meme-software/4.11.0/meme_4.11.0.tar.gz
tar -zxvf meme_4.11.0.tar.gz 1>/dev/null 2>/dev/null
cd meme_4.11.0
./configure --prefix="$HOME/meme_new" 1>/dev/null 2>/dev/null
make install 1>/dev/null 2>/dev/null
cp "$HOME/meme_new/bin/fasta-shuffle-letters" "$HOME/meme_bin/bin/"

