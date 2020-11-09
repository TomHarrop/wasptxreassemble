
graph: graph.svg

graph.svg: Snakefile
	snakemake \
		--rulegraph \
		--forceall \
		| grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
		| dot -Tsvg \
		> graph.svg
		