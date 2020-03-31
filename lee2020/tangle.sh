#!/bin/sh


rm app/Main.hs
listing-tool-exe --latex doc/notes.tex --template doc/template.hs --output app/Main.hs
