# Directory for the RD files
DOCS_DIR=man

# File extensions
RD_EXTENSION = "Rd"
HTML_EXTENSION = "html"

# List of Rd- and HTML-files
Rd_FILES=$(shell find $(DOCS_DIR) -name '*.Rd')
HTML_FILES=$(patsubst %.Rd,%.html,$(Rd_FILES))

# Path to output directory for HTML
OUTPUT_DIR='html_out'

## Runs all commands
all : clean docs $(HTML_FILES)
	# Moving html files
	@mkdir -p $(OUTPUT_DIR) && mv $(DOCS_DIR)/*.$(HTML_EXTENSION) $(OUTPUT_DIR)

## Cleans the data
clean :
	@rm -rf $(OUTPUT_DIR)/*$(HTML_EXTENSION)
	@rm -rf $(DOCS_DIR)/*$(RD_EXTENSION)

## Converts Rd to HTML
%.html: %.Rd
	@Rscript --quiet -e 'rmarkdown::render("$<")'

## Creates Rd documentation
docs:
	@Rscript --quiet -e 'devtools::document()'

## List all files
list_html:
	@echo $(Rd_FILES)
	@printf "\n"
	@echo $(HTML_FILES)

##############################################################################
# Self Documenting Commands                                                  #
##############################################################################

.DEFAULT_GOAL := show-help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
#   * save line in hold space
#   * purge line
#   * Loop:
#       * append newline + line to hold space
#       * go to next line
#       * if line starts with doc comment, strip comment character off and loop
#   * remove target prerequisites
#   * append hold space (+ newline) to line
#   * replace newline plus comments by `---`
#   * print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
