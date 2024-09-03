/home/warxchampion/Desktop/myAnsys/ansys_inc/v222/fluent/bin/fluent-cleanup.pl warxchampion 33839 CLEANUP_EXITING

LOCALHOST=`hostname -s`
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 46021; else ssh warxchampion kill -9 46021; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 46020; else ssh warxchampion kill -9 46020; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 46019; else ssh warxchampion kill -9 46019; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 44705; else ssh warxchampion kill -9 44705; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 44442; else ssh warxchampion kill -9 44442; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 44193; else ssh warxchampion kill -9 44193; fi

rm -f /home/warxchampion/Desktop/Scrap/pl_0_60X604/cleanup-fluent-warxchampion-44442.sh
