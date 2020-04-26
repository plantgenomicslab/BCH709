#!/bin/bash
git pull
chmod 644 -R episodes/*
git add -A
git commit -m "BCH709"
git push
