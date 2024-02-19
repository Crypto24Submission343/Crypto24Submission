fetch_fplll:
	git clone https://github.com/fplll/fplll
	cd fplll; git checkout 5.4.2

patch_fplll:
	patch fplll/fplll/enum/enumerate_base.h diff/fplll_enum_enumerate_base.h.patch
	patch fplll/fplll/enum/enumerate_base.cpp diff/fplll_enum_enumerate_base.cpp.patch
	patch fplll/fplll/enum/enumerate.cpp diff/fplll_enum_enumerate.cpp.patch
	patch fplll/fplll/enum/enumerate.h diff/fplll_enum_enumerate.h.patch

	patch fplll/fplll/pruner/pruner_cost.cpp diff/fplll_pruner_pruner_cost.cpp.patch
	patch fplll/fplll/pruner/pruner_prob.cpp diff/fplll_pruner_pruner_prob.cpp.patch
	patch fplll/fplll/pruner/pruner.cpp diff/fplll_pruner_pruner.cpp.patch
	patch fplll/fplll/pruner/pruner.h diff/fplll_pruner_pruner.h.patch

	patch fplll/tests/test_bkz.cpp diff/tests_test_bkz.cpp.patch
	patch fplll/tests/test_svp.cpp diff/tests_test_svp.cpp.patch
	patch fplll/tests/test_pruner.cpp diff/tests_test_pruner.cpp.patch

	patch fplll/fplll/bkz_param.h diff/fplll_bkz_param.h.patch
	patch fplll/fplll/defs.h diff/fplll_defs.h.patch
	patch fplll/fplll/main.h diff/fplll_main.h.patch
	patch fplll/fplll/main.cpp diff/fplll_main.cpp.patch
	patch fplll/fplll/bkz.h diff/fplll_bkz.h.patch
	patch fplll/fplll/bkz.cpp diff/fplll_bkz.cpp.patch
	patch fplll/fplll/svpcvp.cpp diff/fplll_svpcvp.cpp.patch

	cp diff/fplll_test_gsaalpha.const fplll/tests/gsaalpha.const
	cp diff/fplll_test_gsaalpha.const fplll/fplll/pruner/gsaalpha.const

OUT_PATH=$(PWD)/out
first_compilation:
	-rm -rf ${OUT_PATH}
	mkdir ${OUT_PATH}
	cd fplll; ./autogen.sh;
	cd fplll; ./configure --prefix=${OUT_PATH} --with-max-parallel-enum-dim=0 --with-qd=/usr/lib/x86_64-linux-gnu/
	cd fplll; make -j8
	cd fplll; make install

compile_after_changes:
	cd fplll; make -j8
	cd fplll; make install

FACE_pruning:
	wget https://www2.nict.go.jp/security/pbkzcode/pbkzlib-202205.zip
	unzip pbkzlib-202205.zip -d pbkz
	cd pbkz; bash install_files.sh
	patch pbkz/main.cpp diff/pbkz.patch
	patch pbkz/bkz.conf diff/pbkz_bkz_conf.patch
	patch pbkz/external/tools.h diff/pbkz_external_tools.patch
	cp diff/fplll_test_gsaalpha.const pbkz/gsaalpha.const
	cd pbkz; make all
	cp pbkz/a.out pbkz/FACE_pruning

prepare_experiments:
	-mkdir plots
	-mkdir subtree-plots
	-mkdir stats
	sage se93_strat.py

run_subtree_experiments:
	cat subtree_experiments | parallel --compress

plot_subtree_experiments:
	sed -e 's/$$/ :/' subtree_experiments | parallel --ungroup

run_jensen_experiments:
	cat jensen_experiments | parallel --ungroup

plot_jensen_experiments:
	sed -e 's/$$/ :/' jensen_experiments | parallel --ungroup

generate_patch:
	cd fplll; bash ../gen_patch.sh
	-rm -rf diff
	mkdir diff
	mv fplll/*.patch diff/
	cp fplll/tests/gsaalpha.const diff/fplll_test_gsaalpha.const
	rm diff/configure.ac.patch

clean_bin:
	-rm -rf ${OUT_PATH}

clean_data:
	-rm -rf results
	-rm -rf plots
	-rm -rf subtree-plots
	-rm -rf stats

deep_clean: clean_bin clean_data
	-rm -rf fplll
	-rm -rf pbkz
	-rm -rf no-pruning.json
