#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ezxml.h"
#include "fortprintf.h"

int parse_filters_from_registry(ezxml_t registry);

int parse_filters_from_registry(ezxml_t registry){/*{{{*/
	ezxml_t filters, filter;
	ezxml_t filters2, filter2;
	ezxml_t filter_list, filter_list2;

	const char *filtername, *filtermodule;
	const char *filter2name, *filter2module;
	const char *filterlistname, *filterlistname2;
	const char *corename;

	FILE *fd;

	char core_string[1024];

	int write_use, write_list, wrote_lists;

	corename = ezxml_attr(registry, "core_abbrev");

	sprintf(core_string, "%s", corename);

	fd = fopen("filter_function.inc", "w+");

	fortprintf(fd, "   function %s_setup_filters(inList, filterLists) result(iErr)\n", core_string);

	fortprintf(fd, "      use mpas_filter_list\n");
	fortprintf(fd, "      use mpas_pool_routines\n");

	// Generate use statements for filter setup function
	for (filters = ezxml_child(registry, "filters"); filters; filters = filters->next){
		for (filter = ezxml_child(filters, "filter"); filter; filter = filter->next){
			filtername = ezxml_attr(filter, "name");
			filtermodule = ezxml_attr(filter, "module_name");

			write_use = 1;

			// Search the rest of the current filters block
			for (filter2 = filter->next; filter2; filter2 = filter2->next){
				filter2name = ezxml_attr(filter2, "name");
				filter2module = ezxml_attr(filter2, "module_name");

				if ( strcmp(filtermodule, filter2module) == 0 ) {
					write_use = 0;
				}
			}

			// Search the rest of the filters blocks
			for (filters2 = filters->next; filters2; filters2 = filters2->next){
				for (filter2 = ezxml_child(filters2, "filter"); filter2; filter2 = filter2->next){
					filter2name = ezxml_attr(filter2, "name");
					filter2module = ezxml_attr(filter2, "module_name");

					if ( strcmp(filtermodule, filter2module) == 0 ) {
						write_use = 0;
					}
				}
			}

			if ( write_use ) {
				fortprintf(fd, "      use %s\n", filtermodule);
			}
		}
	}

	fortprintf(fd, "\n");
	fortprintf(fd, "      type (mpas_filter_list_type), pointer :: inList\n");
	fortprintf(fd, "      type (mpas_pool_type), pointer :: filterLists\n");
	fortprintf(fd, "      integer :: iErr\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (mpas_pool_type), pointer :: filterList\n");
	fortprintf(fd, "      procedure (mpas_filter_list_init_function), pointer :: initFunc\n");
	fortprintf(fd, "      procedure (mpas_filter_list_compute_function), pointer :: computeFunc\n");
	fortprintf(fd, "      procedure (mpas_filter_list_restart_function), pointer :: restartFunc\n");
	fortprintf(fd, "      procedure (mpas_filter_list_finalize_function), pointer :: finalizeFunc\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      iErr = 0\n");
	fortprintf(fd, "\n");

	// Register filters properly
	for (filters = ezxml_child(registry, "filters"); filters; filters = filters->next){
		for (filter = ezxml_child(filters, "filter"); filter; filter = filter->next){
			filtername = ezxml_attr(filter, "name");
			filtermodule = ezxml_attr(filter, "module_name");

			fortprintf(fd, "      initFunc => %s_init\n", filtermodule);
			fortprintf(fd, "      computeFunc => %s_compute\n", filtermodule);
			fortprintf(fd, "      restartFunc => %s_restart\n", filtermodule);
			fortprintf(fd, "      finalizeFunc => %s_finalize\n", filtermodule);
			fortprintf(fd, "      call mpas_filter_list_register_filter(inList, '%s', initFunc, computeFunc, restartFunc, finalizeFunc)\n", filtername);
			fortprintf(fd, "\n");

		}
	}

	// Extract lists of filters for proper ordering
	wrote_lists = 0;
	for (filter_list = ezxml_child(registry, "filter_list"); filter_list; filter_list = filter_list->next){
		filterlistname = ezxml_attr(filter_list, "name");
		write_list = 1;
		if ( filter_list->next != NULL ) {
			for ( filter_list2 = filter_list->next; filter_list2; filter_list2 = filter_list2->next){
				filterlistname2 = ezxml_attr(filter_list2, "name");
				if ( strcmp(filterlistname, filterlistname2) == 0 ) {
					write_list = 0;
				}
			}
		}

		if ( write_list ) {
			wrote_lists = 1;
			fortprintf(fd, "\n");
			fortprintf(fd, "      call mpas_pool_create_pool(filterList)\n");

			// Make sure to get filters from all filter_lists with the same name.
			for ( filter_list2 = ezxml_child(registry, "filter_list"); filter_list2; filter_list2 = filter_list2->next){
				filterlistname2 = ezxml_attr(filter_list2, "name");
				if ( strcmp(filterlistname, filterlistname2) == 0 ){
					filter = ezxml_child(filter_list2, "filter");

					for ( filter = ezxml_child(filter_list2, "filter"); filter; filter = filter->next){
						filtername = ezxml_attr(filter, "name");
						fortprintf(fd, "      call mpas_pool_add_config(filterList, '%s', 1)\n", filtername);
					}
				}
			}

			fortprintf(fd, "      call mpas_pool_add_subpool(filterLists, '%s', filterList)\n", filterlistname);
			fortprintf(fd, "\n");
		}
	}

	// If no lists were generated, generate a default list containing all filters in order.
	if ( ! wrote_lists ) {
		fortprintf(fd, "\n");
		fortprintf(fd, "      call mpas_pool_create_pool(filterList)\n");
		for (filter_list = ezxml_child(registry, "filters"); filter_list; filter_list = filter_list->next){
			for (filter = ezxml_child(filter_list, "filter"); filter; filter = filter->next){
				filtername = ezxml_attr(filter, "name");
				fortprintf(fd, "      call mpas_pool_add_config(filterList, '%s', 1)\n", filtername);
			}
		}

		fortprintf(fd, "      call mpas_pool_add_subpool(filterLists, 'default', filterList)\n");
		fortprintf(fd, "\n");
	}

	fortprintf(fd, "   end function %s_setup_filters\n", core_string);


	fclose(fd);

	return 0;
}/*}}}*/

int main(int argc, char ** argv)/*{{{*/
{
	FILE * regfile;
	int err;

	if (argc != 2) {
		fprintf(stderr,"Reading registry file from standard input\n");
		regfile = stdin;
	}
	else if (!(regfile = fopen(argv[1], "r"))) {
		fprintf(stderr,"\nError: Could not open file %s for reading.\n\n", argv[1]);
		return 1;
	}   

	ezxml_t registry = ezxml_parse_fp(regfile);

	if ( parse_filters_from_registry(registry) ) {
		fprintf(stderr, "Parsing failed.....\n");
		return 1;
	}

	return 0;
}/*}}}*/
