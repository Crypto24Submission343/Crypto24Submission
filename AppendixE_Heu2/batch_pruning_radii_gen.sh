# This file runs pruning radii computation for the Kyber-relevant dimensions and single-enumeration probability 2^-64.
# NOTE: All processes will run in parallel.
#
# Before running, `test_pruner` needs to be compiled:
# cd fplll/tests; make test_pruner

bash -c "./fplll/tests/test_pruner -d 406 -p 300 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 0 > d406_s0.json" &
bash -c "./fplll/tests/test_pruner -d 406 -p 300 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 1 > d406_s1.json" &
bash -c "./fplll/tests/test_pruner -d 406 -p 300 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 2 > d406_s2.json" &
bash -c "./fplll/tests/test_pruner -d 406 -p 300 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 3 > d406_s3.json" &

bash -c "./fplll/tests/test_pruner -d 622 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 0 > d622_s0.json" &
bash -c "./fplll/tests/test_pruner -d 622 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 1 > d622_s1.json" &
bash -c "./fplll/tests/test_pruner -d 622 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 2 > d622_s2.json" &
bash -c "./fplll/tests/test_pruner -d 622 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 3 > d622_s3.json" &

bash -c "./fplll/tests/test_pruner -d 623 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 0 > d623_s0.json" &
bash -c "./fplll/tests/test_pruner -d 623 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 1 > d623_s1.json" &
bash -c "./fplll/tests/test_pruner -d 623 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 2 > d623_s2.json" &
bash -c "./fplll/tests/test_pruner -d 623 -p 350 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 3 > d623_s3.json" &

bash -c "./fplll/tests/test_pruner -d 872 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 0 > d872_s0.json" &
bash -c "./fplll/tests/test_pruner -d 872 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 1 > d872_s1.json" &
bash -c "./fplll/tests/test_pruner -d 872 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 2 > d872_s2.json" &
bash -c "./fplll/tests/test_pruner -d 872 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 3 > d872_s3.json" &

bash -c "./fplll/tests/test_pruner -d 873 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 0 > d873_s0.json" &
bash -c "./fplll/tests/test_pruner -d 873 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 1 > d873_s1.json" &
bash -c "./fplll/tests/test_pruner -d 873 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 2 > d873_s2.json" &
bash -c "./fplll/tests/test_pruner -d 873 -p 500 -t 5.42101086242752217003726400434970855712890625e-20 -j -s 3 > d873_s3.json" &
