all: clean install-report
	@echo "Compile the quartet-dseqc-report...."
	@bin/lein uberjar
	@printf "\n\n\e[1;32mRun the command for more details: \nsource .env/bin/activate\njava -jar target/uberjar/quartet-dseqc-report-*-standalone.jar -h\e[0m"

clean:
	@echo "Clean the environment..."
	@bin/lein clean
	@rm -rf .env .lsp .clj-kondo report/dist report/quartet-proteome-report.egg-info exp2qcdt.tar.gz resources/renv/library resources/renv/staging

make-env:
	virtualenv -p python3 .env

install-report: make-env
	. .env/bin/activate
	cd report && python3 setup.py sdist && pip3 install dist/*.tar.gz
