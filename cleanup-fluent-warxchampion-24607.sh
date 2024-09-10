/home/warxchampion/Desktop/myAnsys/ansys_inc/v222/fluent/bin/fluent-cleanup.pl warxchampion 33353 CLEANUP_EXITING

LOCALHOST=`hostname -s`
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24928; else ssh warxchampion kill -9 24928; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24927; else ssh warxchampion kill -9 24927; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24926; else ssh warxchampion kill -9 24926; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24925; else ssh warxchampion kill -9 24925; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24924; else ssh warxchampion kill -9 24924; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24923; else ssh warxchampion kill -9 24923; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24922; else ssh warxchampion kill -9 24922; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24921; else ssh warxchampion kill -9 24921; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24920; else ssh warxchampion kill -9 24920; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24919; else ssh warxchampion kill -9 24919; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24607; else ssh warxchampion kill -9 24607; fi
if [[ warxchampion == "$LOCALHOST"* ]]; then kill -9 24379; else ssh warxchampion kill -9 24379; fi

rm -f /home/warxchampion/Desktop/FoamDATA/GitData/FTELabWork/cleanup-fluent-warxchampion-24607.sh
