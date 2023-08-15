virtual_env:
	poetry env use python3.6
	poetry shell

update_deps:
	poetry install
	poetry lock

.PHONY: usepy36
