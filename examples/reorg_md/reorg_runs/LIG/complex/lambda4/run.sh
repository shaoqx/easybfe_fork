#!/usr/bin/env bash

cd 01.em
echo Running 01.em && touch running.tag
if [ ! -f done.tag ]; then
  source 01.em.sh > 01.em.stdout 2>&1
  if [ $? -ne 0 ]; then
    mv running.tag error.tag && echo "Error occurs!" && exit 1
  fi
  mv running.tag done.tag
fi
cd ..

cd 02.heat
echo Running 02.heat && touch running.tag
if [ ! -f done.tag ]; then
  source 02.heat.sh > 02.heat.stdout 2>&1
  if [ $? -ne 0 ]; then
    mv running.tag error.tag && echo "Error occurs!" && exit 1
  fi
  mv running.tag done.tag
fi
cd ..

cd 03.pres
echo Running 03.pres && touch running.tag
if [ ! -f done.tag ]; then
  source 03.pres.sh > 03.pres.stdout 2>&1
  if [ $? -ne 0 ]; then
    mv running.tag error.tag && echo "Error occurs!" && exit 1
  fi
  mv running.tag done.tag
fi
cd ..

cd 04.pre_prod
echo Running 04.pre_prod && touch running.tag
if [ ! -f done.tag ]; then
  source 04.pre_prod.sh > 04.pre_prod.stdout 2>&1
  if [ $? -ne 0 ]; then
    mv running.tag error.tag && echo "Error occurs!" && exit 1
  fi
  mv running.tag done.tag
fi
cd ..

cd 05.prod
echo Running 05.prod && touch running.tag
if [ ! -f done.tag ]; then
  source 05.prod.sh > 05.prod.stdout 2>&1
  if [ $? -ne 0 ]; then
    mv running.tag error.tag && echo "Error occurs!" && exit 1
  fi
  mv running.tag done.tag
fi
cd ..

