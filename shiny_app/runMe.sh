#!/bin/bash

port="8083"

for i in $(ps aux | grep 'shiny::runApp' | grep ${port} | awk '{print $2}'); do
  kill -9 ${i}
done

git pull &&
  nohup R -e "shiny::runApp('.', port = ${port}, host = '0.0.0.0')" \
    >runMe.log &
