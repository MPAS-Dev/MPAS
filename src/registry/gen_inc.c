// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ezxml/ezxml.h"
#include "registry_types.h"
#include "gen_inc.h"
#include "fortprintf.h"
#include "utility.h"

#define STR(s) #s
#define MACRO_TO_STR(s) STR(s)

void write_model_variables(ezxml_t registry){/*{{{*/
	const char * suffix = MACRO_TO_STR(MPAS_NAMELIST_SUFFIX);
	const char * exe_name = MACRO_TO_STR(MPAS_EXE_NAME);
	const char * git_ver = MACRO_TO_STR(MPAS_GIT_VERSION);

	const char *modelname, *corename, *version;
	FILE *fd;

	modelname = ezxml_attr(registry, "model");
	corename = ezxml_attr(registry, "core");
	version = ezxml_attr(registry, "version");

	fd = fopen("model_variables.inc", "w+");

	fortprintf(fd, "       character (len=StrKIND) :: modelName = '%s' !< Constant: Name of model\n", modelname);
	fortprintf(fd, "       character (len=StrKIND) :: coreName = '%s' !< Constant: Name of core\n", corename);
	fortprintf(fd, "       character (len=StrKIND) :: modelVersion = '%s' !< Constant: Version number\n", version);
	fortprintf(fd, "       character (len=StrKIND) :: namelist_filename = 'namelist.%s' !< Constant: Name of namelist file\n", suffix);
	fortprintf(fd, "       character (len=StrKIND) :: streams_filename = 'streams.%s' !< Constant: Name of stream configuration file\n", suffix);
	fortprintf(fd, "       character (len=StrKIND) :: executableName = '%s' !< Constant: Name of executable generated at build time.\n", exe_name);
	fortprintf(fd, "       character (len=StrKIND) :: git_version = '%s' !< Constant: Version string from git-describe.\n", git_ver);

	fclose(fd);
}/*}}}*/


int write_field_pointers(FILE* fd){/*{{{*/
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (field0DReal), pointer :: r0Ptr\n");
	fortprintf(fd, "      type (field1DReal), pointer :: r1Ptr\n");
	fortprintf(fd, "      type (field2DReal), pointer :: r2Ptr\n");
	fortprintf(fd, "      type (field3DReal), pointer :: r3Ptr\n");
	fortprintf(fd, "      type (field4DReal), pointer :: r4Ptr\n");
	fortprintf(fd, "      type (field5DReal), pointer :: r5Ptr\n");
	fortprintf(fd, "      type (field0DInteger), pointer :: i0Ptr\n");
	fortprintf(fd, "      type (field1DInteger), pointer :: i1Ptr\n");
	fortprintf(fd, "      type (field2DInteger), pointer :: i2Ptr\n");
	fortprintf(fd, "      type (field3DInteger), pointer :: i3Ptr\n");
	fortprintf(fd, "      type (field0DChar), pointer :: c0Ptr\n");
	fortprintf(fd, "      type (field1DChar), pointer :: c1Ptr\n");
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int write_field_pointer_arrays(FILE* fd){/*{{{*/
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (field0DReal), dimension(:), pointer :: r0Ptr\n");
	fortprintf(fd, "      type (field1DReal), dimension(:), pointer :: r1Ptr\n");
	fortprintf(fd, "      type (field2DReal), dimension(:), pointer :: r2Ptr\n");
	fortprintf(fd, "      type (field3DReal), dimension(:), pointer :: r3Ptr\n");
	fortprintf(fd, "      type (field4DReal), dimension(:), pointer :: r4Ptr\n");
	fortprintf(fd, "      type (field5DReal), dimension(:), pointer :: r5Ptr\n");
	fortprintf(fd, "      type (field0DInteger), dimension(:), pointer :: i0Ptr\n");
	fortprintf(fd, "      type (field1DInteger), dimension(:), pointer :: i1Ptr\n");
	fortprintf(fd, "      type (field2DInteger), dimension(:), pointer :: i2Ptr\n");
	fortprintf(fd, "      type (field3DInteger), dimension(:), pointer :: i3Ptr\n");
	fortprintf(fd, "      type (field0DChar), dimension(:), pointer :: c0Ptr\n");
	fortprintf(fd, "      type (field1DChar), dimension(:), pointer :: c1Ptr\n");
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int set_pointer_name(int type, int ndims, char *pointer_name){/*{{{*/
	if(type == REAL) {
		switch (ndims){
			default:
			case 0:
				snprintf(pointer_name, 1024, "r0Ptr");
				break;
			case 1:
				snprintf(pointer_name, 1024, "r1Ptr");
				break;
			case 2:
				snprintf(pointer_name, 1024, "r2Ptr");
				break;
			case 3:
				snprintf(pointer_name, 1024, "r3Ptr");
				break;
			case 4:
				snprintf(pointer_name, 1024, "r4Ptr");
				break;
			case 5:
				snprintf(pointer_name, 1024, "r5Ptr");
				break;
		}
	} else if (type == INTEGER) {
		switch (ndims){
			default:
			case 0:
				snprintf(pointer_name, 1024, "i0Ptr");
				break;
			case 1:
				snprintf(pointer_name, 1024, "i1Ptr");
				break;
			case 2:
				snprintf(pointer_name, 1024, "i2Ptr");
				break;
			case 3:
				snprintf(pointer_name, 1024, "i3Ptr");
				break;
		}
	} else if (type == CHARACTER) {
		switch (ndims){
			default:
			case 0:
				snprintf(pointer_name, 1024, "c0Ptr");
				break;
			case 1:
				snprintf(pointer_name, 1024, "c1Ptr");
				break;
		}
	}

	return 0;
}/*}}}*/


int add_package_to_list(const char * package, const char * package_list){/*{{{*/
	char *token, *string, *tofree;

	string = strdup(package_list);
	tofree = string;
	token = strsep(&string, ";");

	if(strcmp(package, token) == 0){
		return 0;
	}

	while( (token = strsep(&string, ";")) != NULL){
		if(strcmp(package, token) == 0){

			return 0;
		}
	}

	return 1;
}/*}}}*/


int build_struct_package_lists(ezxml_t currentPosition, char * out_packages){/*{{{*/
	ezxml_t child_xml1, child_xml2;

	const char *package_list;
	const char *name;

	char *token, *string, *tofree;
	int empty_packages;
	int empty_struct;

	package_list = ezxml_attr(currentPosition, "packages");

	empty_packages = 0;
	empty_struct = 1;

	// Check for vars that don't have packages.
	for(child_xml1 = ezxml_child(currentPosition, "var"); child_xml1 && !empty_packages; child_xml1 = child_xml1->next){
		package_list = ezxml_attr(child_xml1, "packages");

		if(!package_list){
			empty_packages = 1;
		}
		empty_struct = 0;
	}

	// Check for vararrays and constituents that don't have packages.
	for(child_xml1 = ezxml_child(currentPosition, "var_array"); child_xml1 && !empty_packages; child_xml1 = child_xml1->next){
		package_list = ezxml_attr(child_xml1, "packages");

		if(!package_list){
			for(child_xml2 = ezxml_child(child_xml1, "var"); child_xml2 && !empty_packages; child_xml2 = child_xml2->next){
				package_list = ezxml_attr(child_xml2, "packages");

				if(!package_list){
					empty_packages = 1;
				}
				empty_struct = 0;
			}
		}
	}

	// If any var/var_array doesn't have packages on it, the struct doesn't have packages on it.
	if(empty_packages || empty_struct){
		return 1;
	} else {
		// Build unique list of packages from nested vars and var arrays.
		for(child_xml1 = ezxml_child(currentPosition, "var_array"); child_xml1; child_xml1 = child_xml1->next){
			package_list = ezxml_attr(child_xml1, "packages");

			// Build list of unique packages from var_array
			if(package_list){
				string = strdup(package_list);
				tofree = string;
				token = strsep(&string, ";");

				if(out_packages[0] == '\0'){
					sprintf(out_packages, "%s", token);
				} else if(add_package_to_list(token, out_packages)){
					sprintf(out_packages, "%s;%s", out_packages, token);
				}

				while( (token = strsep(&string, ";")) != NULL){
					if(add_package_to_list(token, out_packages)){
						sprintf(out_packages, "%s;%s", out_packages, token);
					}
				}

				free(tofree);
			}

			for(child_xml2 = ezxml_child(child_xml1, "var"); child_xml2; child_xml2 = child_xml2->next){
				package_list = ezxml_attr(child_xml2, "packages");

				// Build list of unique packages from child var
				if(package_list){
					string = strdup(package_list);
					tofree = string;
					token = strsep(&string, ";");

					if(out_packages[0] == '\0'){
						sprintf(out_packages, "%s", token);
					} else if(add_package_to_list(token, out_packages)){
						sprintf(out_packages, "%s;%s", out_packages, token);
					}

					while( (token = strsep(&string, ";")) != NULL){
						if(add_package_to_list(token, out_packages)){
							sprintf(out_packages, "%s;%s", out_packages, token);
						}
					}

					free(tofree);
				}
			}
		}

		for(child_xml1 = ezxml_child(currentPosition, "var"); child_xml1; child_xml1 = child_xml1->next){
			package_list = ezxml_attr(child_xml1, "packages");

			// Build list of unique packages from child var
			if(package_list){
				string = strdup(package_list);
				tofree = string;
				token = strsep(&string, ";");

				if(out_packages[0] == '\0'){
					sprintf(out_packages, "%s", token);
				} else if(add_package_to_list(token, out_packages)){
					sprintf(out_packages, "%s;%s", out_packages, token);
				}

				while( (token = strsep(&string, ";")) != NULL){
					if(add_package_to_list(token, out_packages)){
						sprintf(out_packages, "%s;%s", out_packages, token);
					}
				}

				free(tofree);
			}
		}
		return 0;
	}
}/*}}}*/


int get_dimension_information(const char *dims, int *ndims, int *has_time, int *decomp){/*{{{*/
	char *string, *tofree, *token;

	(*ndims) = 0;
	(*decomp) = -1;
	(*has_time) = 0;
	string = strdup(dims);
	tofree = string;
	token = strsep(&string, " ");
	if(strlen(token) > 0){
		if(strcmp(token, "Time") != 0){
			if(strcmp(token, "nCells") == 0){
				if((*decomp) == -1){
					(*decomp) = CELLS;
				} else {
					printf("ERROR: Multiple decomposition types\n");
					return 1;
				}
			} else if (strcmp(token, "nEdges") == 0){
				if((*decomp) == -1){
					(*decomp) = EDGES;
				} else {
					printf("ERROR: Multiple decomposition types\n");
					return 1;
				}
			} else if (strcmp(token, "nVertices") == 0){
				if((*decomp) == -1){
					(*decomp) = VERTICES;
				} else {
					printf("ERROR: Multiple decomposition types\n");
					return 1;
				}
			}
			(*ndims)++;
		} else {
			(*has_time) = 1;
		}
	}
	while( (token = strsep(&string, " ")) != NULL){
		if(strcmp(token, "Time") != 0){
			if(strcmp(token, "nCells") == 0){
				if((*decomp) == -1){
					(*decomp) = CELLS;
				} else {
					printf("ERROR: Multiple decomposition types\n");
					return 1;
				}
			} else if (strcmp(token, "nEdges") == 0){
				if((*decomp) == -1){
					(*decomp) = EDGES;
				} else {
					printf("ERROR: Multiple decomposition types\n");
					return 1;
				}
			} else if (strcmp(token, "nVertices") == 0){
				if((*decomp) == -1){
					(*decomp) = VERTICES;
				} else {
					printf("ERROR: Multiple decomposition types\n");
					return 1;
				}
			}
			(*ndims)++;
		} else {
			(*has_time) = 1;
		}
	}
	free(tofree);

	return 0;
}/*}}}*/


int get_field_information(const char *vartype, const char *varval, char *default_value, int *type){/*{{{*/
	if (strcmp(vartype, "real") == 0){
		(*type) = REAL;
		if(!varval){
			snprintf(default_value, 1024, "0.0");
		} else {
			snprintf(default_value, 1024, "%s", varval);
		}
	} else if (strcmp(vartype, "integer") == 0){
		(*type) = INTEGER;
		if(!varval){
			snprintf(default_value, 1024, "0");
		} else {
			snprintf(default_value, 1024, "%s", varval);
		}
	} else if (strcmp(vartype, "text") == 0){
		(*type) = CHARACTER;
		if(!varval){
			snprintf(default_value, 1024, "''");
		} else {
			snprintf(default_value, 1024, "%s", varval);
		}
	}

	return 0;
}/*}}}*/


void write_default_namelist(ezxml_t registry) /*{{{*/
{
	ezxml_t rec_xml, opt_xml;

	const char *recname, *optname, *opttype, *optdefault;
	const char *recindef, *optindef;
	FILE *fd;
	char filename[1024];
	const char * suffix = MACRO_TO_STR(MPAS_NAMELIST_SUFFIX);
	int print_record, print_option;


	sprintf(filename, "namelist.%s.defaults", suffix);
	fd = fopen(filename, "w+");

	for(rec_xml = ezxml_child(registry, "nml_record"); rec_xml; rec_xml = rec_xml->next){
		recname = ezxml_attr(rec_xml, "name");
		recindef = ezxml_attr(rec_xml, "in_defaults");

		print_record = 1;

		if(recindef){
			if(strcmp(recindef, "true") != 0){
				print_record = 0;
			}
		}

		if(print_record) {
			fprintf(fd, "&%s\n", recname);
			for(opt_xml = ezxml_child(rec_xml, "nml_option"); opt_xml; opt_xml = opt_xml->next){
				optname = ezxml_attr(opt_xml, "name");
				opttype = ezxml_attr(opt_xml, "type");
				optdefault = ezxml_attr(opt_xml, "default_value");
				optindef = ezxml_attr(opt_xml, "in_defaults");

				print_option = 1;

				if (optindef) {
					if (strcmp(optindef, "true") != 0){
						print_option = 0;
					}
				}

				if(print_option){
					if (strcmp(opttype, "real") == 0){
						fprintf(fd, "    %s = %s\n", optname, optdefault);
					} else if (strcmp(opttype, "integer") == 0){
						fprintf(fd, "    %s = %s\n", optname, optdefault);
					} else if (strcmp(opttype, "logical") == 0){
						if (strcmp(optdefault, "true") == 0 || strcmp(optdefault, ".true.") == 0){
							fprintf(fd, "    %s = .true.\n", optname);
						} else {
							fprintf(fd, "    %s = .false.\n", optname);
						}
					} else if (strcmp(opttype, "character") == 0){
						fprintf(fd, "    %s = '%s'\n", optname, optdefault);
					}
				}
			}
			fprintf(fd, "/\n\n");
		}
	}

	fclose(fd);
}/*}}}*/


void write_default_streams(ezxml_t registry){/*{{{*/
	ezxml_t streams_xml, varstruct_xml, opt_xml, var_xml, vararray_xml, stream_xml;

	const char *optstream, *optname, *optvarname, *opttype;
	const char *optimmutable, *optfilename, *optfilename_interval, *optinterval_in, *optinterval_out, *optruntime, *optpackages;
	const char *optref_time, *optprecision, *optrecord_interval, *optclobber_mode;
	FILE *fd, *fd2;
	char filename[64], filename2[64];
	const char * suffix = MACRO_TO_STR(MPAS_NAMELIST_SUFFIX);


	sprintf(filename, "streams.%s.defaults", suffix);
	fd = fopen(filename, "w");

	fprintf(stderr, "Generating run-time stream definitions\n");

	fprintf(fd, "<streams>\n\n");
	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (opt_xml = ezxml_child(streams_xml, "stream"); opt_xml; opt_xml = opt_xml->next) {
			optname = ezxml_attr(opt_xml, "name");
			opttype = ezxml_attr(opt_xml, "type");
			optfilename = ezxml_attr(opt_xml, "filename_template");
			optfilename_interval = ezxml_attr(opt_xml, "filename_interval");
			optinterval_in = ezxml_attr(opt_xml, "input_interval");
			optinterval_out = ezxml_attr(opt_xml, "output_interval");
			optimmutable = ezxml_attr(opt_xml, "immutable");
			optruntime = ezxml_attr(opt_xml, "runtime_format");
			optpackages = ezxml_attr(opt_xml, "packages");
			optprecision = ezxml_attr(opt_xml, "precision");
			optref_time = ezxml_attr(opt_xml, "reference_time");
			optrecord_interval = ezxml_attr(opt_xml, "record_interval");
			optclobber_mode = ezxml_attr(opt_xml, "clobber_mode");

			/* Generate immutable default stream */
			if (optimmutable != NULL && strcmp(optimmutable, "true") == 0) {
				fprintf(fd, "<immutable_stream name=\"%s\"\n", optname);
				fprintf(fd, "                  type=\"%s\"\n", opttype);
				fprintf(fd, "                  filename_template=\"%s\"\n", optfilename);
				if (optfilename_interval) {
					fprintf(fd, "                  filename_interval=\"%s\"\n", optfilename_interval);
				}
				if (optpackages) {
					fprintf(fd, "                  packages=\"%s\"\n",optpackages);
				}
				if (optprecision) {
					fprintf(fd, "                  precision=\"%s\"\n", optprecision);
				}
				if (optref_time) {
					fprintf(fd, "                  reference_time=\"%s\"\n", optref_time);
				}
				if (optrecord_interval) {
					fprintf(fd, "                  record_interval=\"%s\"\n", optrecord_interval);
				}
				if (optclobber_mode) {
					fprintf(fd, "                  clobber_mode=\"%s\"\n", optclobber_mode);
				}
				if (strstr(opttype, "input") != NULL) {
					fprintf(fd, "                  input_interval=\"%s\"", optinterval_in);
				}
				if (strstr(opttype, "output") != NULL) {
					if (strstr(opttype, "input") != NULL) fprintf(fd, "\n");
					fprintf(fd, "                  output_interval=\"%s\"", optinterval_out);
				}
				fprintf(fd, "/>\n\n");
			}
			else {
				fprintf(fd, "<stream name=\"%s\"\n", optname);
				fprintf(fd, "        type=\"%s\"\n", opttype);
				fprintf(fd, "        filename_template=\"%s\"\n", optfilename);
				if( optfilename_interval) {
					fprintf(fd, "        filename_interval=\"%s\"\n", optfilename_interval);
				}
				if (optpackages) {
					fprintf(fd, "        packages=\"%s\"\n",optpackages);
				}
				if (optprecision) {
					fprintf(fd, "        precision=\"%s\"\n", optprecision);
				}
				if (optref_time) {
					fprintf(fd, "        reference_time=\"%s\"\n", optref_time);
				}
				if (optrecord_interval) {
					fprintf(fd, "        record_interval=\"%s\"\n", optrecord_interval);
				}
				if (optclobber_mode) {
					fprintf(fd, "        clobber_mode=\"%s\"\n", optclobber_mode);
				}
				if (strstr(opttype, "input") != NULL) {
					fprintf(fd, "        input_interval=\"%s\"", optinterval_in);
				}
				if (strstr(opttype, "output") != NULL) {
					if (strstr(opttype, "input") != NULL) fprintf(fd, "\n");
					fprintf(fd, "        output_interval=\"%s\"", optinterval_out);
				}
				fprintf(fd, ">\n\n");

				/*
				 * Depending on the runtime format, we either generate a separate list of fields for
				 *   each stream, or we list the fields directly in the main stream control file
				 */

				if (strcmp(optruntime,"single_file") == 0) {

					/* Loop over streams listed within the stream */
					for (stream_xml = ezxml_child(opt_xml, "stream"); stream_xml; stream_xml = stream_xml->next){
						optname = ezxml_attr(stream_xml, "name");
						fprintf(fd, "    <stream name=\"%s\"/>\n", optname);
					}

					/* Loop over fields looking for any that belong to the stream */
					for (varstruct_xml = ezxml_child(registry, "var_struct"); varstruct_xml; varstruct_xml = varstruct_xml->next) {
						for (var_xml = ezxml_child(varstruct_xml, "var"); var_xml; var_xml = var_xml->next) {
							optstream = ezxml_attr(var_xml, "streams");
							if (optstream != NULL && strstr(optstream, optname) != NULL) {
								optvarname = ezxml_attr(var_xml, "name");
								fprintf(fd, "    <var name=\"%s\"/>\n", optvarname);
							}
						}
					}

					/* Loop over var_arrays listed within the stream */
					for (vararray_xml = ezxml_child(opt_xml, "var_array"); vararray_xml; vararray_xml = vararray_xml->next){
						optname = ezxml_attr(vararray_xml, "name");
						fprintf(fd, "    <var_array name=\"%s\"/>\n", optname);
					}

					/* Loop over fields listed within the stream */
					for (var_xml = ezxml_child(opt_xml, "var"); var_xml; var_xml = var_xml->next) {
						optname = ezxml_attr(var_xml, "name");
						fprintf(fd, "    <var name=\"%s\"/>\n", optname);
					}
				}
				else if (strcmp(optruntime,"separate_file") == 0) {

					sprintf(filename2, "stream_list.%s.%s", suffix, optname);

					fprintf(fd, "    <file name=\"%s\"/>\n", filename2);

					fd2 = fopen(filename2, "w+");

					/* Loop over streams listed within the stream */
					for (stream_xml = ezxml_child(opt_xml, "stream"); stream_xml; stream_xml = stream_xml->next){
						optname = ezxml_attr(stream_xml, "name");
						fprintf(fd, "    <stream name=\"%s\"/>\n", optname);
					}

					/* Loop over fields looking for any that belong to the stream */
					for (varstruct_xml = ezxml_child(registry, "var_struct"); varstruct_xml; varstruct_xml = varstruct_xml->next) {
						for (var_xml = ezxml_child(varstruct_xml, "var"); var_xml; var_xml = var_xml->next) {
							optstream = ezxml_attr(var_xml, "streams");
							if (optstream != NULL && strstr(optstream, optname) != NULL) {
								optvarname = ezxml_attr(var_xml, "name");
								fprintf(fd2, "%s\n", optvarname);
							}
						}
					}

					/* Loop over var_arrays listed within the stream */
					for (vararray_xml = ezxml_child(opt_xml, "var_array"); vararray_xml; vararray_xml = vararray_xml->next){
						optname = ezxml_attr(vararray_xml, "name");
						fprintf(fd2, "%s\n", optname);
					}

					/* Loop over fields listed within the stream */
					for (var_xml = ezxml_child(opt_xml, "var"); var_xml; var_xml = var_xml->next) {
						optname = ezxml_attr(var_xml, "name");
						fprintf(fd2, "%s\n", optname);
					}

					fclose(fd2);

				}
				else {
					fprintf(stderr, "******************************************************\n");
					fprintf(stderr, "Error in specification of stream_format; this probably \n");
					fprintf(stderr, "should have been caught during validation...\n");
					fprintf(stderr, "******************************************************\n");
				}

				fprintf(fd, "\n");
				fprintf(fd, "</stream>\n\n");
			}
		}
	}
	fprintf(fd, "</streams>\n");

	fclose(fd);
}/*}}}*/


int parse_packages_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t packages_xml, package_xml;

	const char *packagename, *packagedesc, *const_core;
	FILE *fd;
	char core_string[1024];

	const_core = ezxml_attr(registry, "core");

	fd = fopen("define_packages.inc", "w+");

	sprintf(core_string, "_%s_", const_core);

	// For now, don't include core name in subroutines
	sprintf(core_string, "_");


	fortprintf(fd, "   subroutine mpas_generate%spackages(packagePool)\n", core_string);
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: packagePool !< Input: MPAS Pool for containing package logicals.\n\n");

	// Parse Packages
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");
			packagedesc = ezxml_attr(package_xml, "description");

			fortprintf(fd, "      call mpas_pool_add_package(packagePool, '%sActive', .false.)\n", packagename);
		}
	}

	fortprintf(fd, "   end subroutine mpas_generate%spackages\n", core_string);

	fclose(fd);

	return 0;
}/*}}}*/


int parse_namelist_records_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t nmlrecs_xml, nmlopt_xml;

	const char *const_core;
	const char *nmlrecname, *nmlrecindef, *nmlrecinsub;
	const char *nmloptname, *nmlopttype, *nmloptval, *nmloptunits, *nmloptdesc, *nmloptposvals, *nmloptindef;

	char pool_name[1024];
	char core_string[1024];

	int in_subpool;

	FILE *fd, *fd2;

	const_core = ezxml_attr(registry, "core");

	sprintf(core_string, "_%s_", const_core);

	// For now, don't include core name in subroutines
	sprintf(core_string, "_");

	fd = fopen("namelist_defines.inc", "w+");
	fd2 = fopen("namelist_call.inc", "w+");

	fortprintf(fd2, "   subroutine mpas_setup_namelists(configPool, namelistFilename, dminfo)\n");
	fortprintf(fd2, "      type (mpas_pool_type), intent(inout) :: configPool\n");
	fortprintf(fd2, "      character (len=*), intent(in) :: namelistFilename\n");
	fortprintf(fd2, "      type (dm_info), intent(in) :: dminfo\n");
	fortprintf(fd2, "\n");
	fortprintf(fd2, "      integer :: unitNumber\n");
	fortprintf(fd2, "\n");
	fortprintf(fd2, "      unitNumber = 21\n");
	fortprintf(fd2, "      open(unitNumber,file=trim(namelistFilename),status='old',form='formatted')\n");
	fortprintf(fd2, "\n");


	// Parse Namelist Records
	for (nmlrecs_xml = ezxml_child(registry, "nml_record"); nmlrecs_xml; nmlrecs_xml = nmlrecs_xml->next){
		nmlrecname = ezxml_attr(nmlrecs_xml, "name");
		nmlrecindef = ezxml_attr(nmlrecs_xml, "in_defaults");
		nmlrecinsub = ezxml_attr(nmlrecs_xml, "in_subpool");

		in_subpool = 0;

		if(nmlrecinsub){
			if(strcmp(nmlrecinsub, "true") == 0){
				in_subpool = 1;
			}
		}

		if(in_subpool){
			sprintf(pool_name, "recordPool");
		} else {
			sprintf(pool_name, "configPool");
		}

		// Add call to driver routine.
		fortprintf(fd2, "      call mpas_setup%snmlrec_%s(configPool, unitNumber, dminfo)\n", core_string, nmlrecname);

		// Start defining new subroutine for namelist record.
		fortprintf(fd, "   subroutine mpas_setup%snmlrec_%s(configPool, unitNumber, dminfo)\n", core_string, nmlrecname);
		fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: configPool\n");
		fortprintf(fd, "      integer, intent(in) :: unitNumber\n");
		fortprintf(fd, "      type (dm_info), intent(in) :: dminfo\n");
		fortprintf(fd, "\n");
		fortprintf(fd, "      integer :: ierr\n");
		fortprintf(fd, "      type (mpas_pool_type), pointer :: recordPool\n");
		fortprintf(fd, "\n");

		// Define variable defintions prior to reading the namelist in.
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");
			nmloptval = ezxml_attr(nmlopt_xml, "default_value");
			nmloptunits = ezxml_attr(nmlopt_xml, "units");
			nmloptdesc = ezxml_attr(nmlopt_xml, "description");
			nmloptposvals = ezxml_attr(nmlopt_xml, "possible_values");
			nmloptindef = ezxml_attr(nmlopt_xml, "in_defaults");

			if(strncmp(nmlopttype, "real", 1024) == 0){
				fortprintf(fd, "      real (kind=RKIND) :: %s = %lf\n", nmloptname, (float)atof(nmloptval));
			} else if(strncmp(nmlopttype, "integer", 1024) == 0){
				fortprintf(fd, "      integer :: %s = %d\n", nmloptname, atoi(nmloptval));
			} else if(strncmp(nmlopttype, "logical", 1024) == 0){
				if(strncmp(nmloptval, "true", 1024) == 0 || strncmp(nmloptval, ".true.", 1024) == 0){
					fortprintf(fd, "      logical :: %s = .true.\n", nmloptname);
				} else {
					fortprintf(fd, "      logical :: %s = .false.\n", nmloptname);
				}
			} else if(strncmp(nmlopttype, "character", 1024) == 0){
					fortprintf(fd, "      character (len=StrKIND) :: %s = '%s'\n", nmloptname, nmloptval);
			}
		}
		fortprintf(fd, "\n");

		// Define the namelist block, to read the namelist record in.
		fortprintf(fd, "      namelist /%s/ &\n", nmlrecname);
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			if(nmlopt_xml->next){
				fortprintf(fd, "         %s, &\n", nmloptname);
			} else {
				fortprintf(fd, "         %s\n", nmloptname);
			}
		}

		if(in_subpool){
			fortprintf(fd, "\n");
			fortprintf(fd, "      allocate(recordPool)\n");
			fortprintf(fd, "      call mpas_pool_create_pool(recordPool)\n");
			fortprintf(fd, "      call mpas_pool_add_subpool(configPool, '%s', recordPool)\n", nmlrecname);
			fortprintf(fd, "\n");
		}

		fortprintf(fd, "      if (dminfo %% my_proc_id == IO_NODE) then\n");
		fortprintf(fd, "         rewind(unitNumber)\n");
		fortprintf(fd, "         read(unitNumber, %s, iostat=ierr)\n", nmlrecname);
		fortprintf(fd, "         if (ierr > 0) then\n");
		fortprintf(fd, "            write(stderrUnit, *) 'Error while reading namelist record %s.'\n", nmlrecname);
		fortprintf(fd, "            call mpas_dmpar_abort(dminfo)\n");
		fortprintf(fd, "         else if (ierr < 0) then\n");
		fortprintf(fd, "            write(stderrUnit,*) 'Namelist record %s not found; using default values for variables in this namelist'\n", nmlrecname);
		fortprintf(fd, "         end if\n");
		fortprintf(fd, "      end if\n");

		// Broadcast ierr, to check if a broadcast should happen for the options (if namelist was read in)
		fortprintf(fd, "      call mpas_dmpar_bcast_int(dminfo, ierr)\n");

		fortprintf(fd, "\n");
		// Define broadcast calls for namelist values.
		fortprintf(fd, "      if (ierr == 0) then\n");
		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");
			nmlopttype = ezxml_attr(nmlopt_xml, "type");

			if(strncmp(nmlopttype, "real", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_real(dminfo, %s)\n", nmloptname);
			} else if(strncmp(nmlopttype, "integer", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_int(dminfo, %s)\n", nmloptname);
			} else if(strncmp(nmlopttype, "logical", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_logical(dminfo, %s)\n", nmloptname);
			} else if(strncmp(nmlopttype, "character", 1024) == 0){
				fortprintf(fd, "         call mpas_dmpar_bcast_char(dminfo, %s)\n", nmloptname);
			}
		}
		fortprintf(fd, "      end if\n");
		fortprintf(fd, "\n");

		for (nmlopt_xml = ezxml_child(nmlrecs_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
			nmloptname = ezxml_attr(nmlopt_xml, "name");

			fortprintf(fd, "      call mpas_pool_add_config(%s, '%s', %s)\n", pool_name, nmloptname, nmloptname);
		}
		fortprintf(fd, "\n");

		// End new subroutine for namelist record.
		fortprintf(fd, "   end subroutine mpas_setup%snmlrec_%s\n", core_string, nmlrecname);
		fortprintf(fd, "\n\n");
	}

	fortprintf(fd2, "\n");
	fortprintf(fd2, "      close(unitNumber)\n");
	fortprintf(fd2, "   end subroutine mpas_setup_namelists\n");

	return 0;
}/*}}}*/


int parse_dimensions_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t nmlrec_xml, nmlopt_xml;

	const char *nmlrecname, *nmlrecinsub, *nmloptname, *nmlopttype;
	const char *dimname, *dimunits, *dimdesc, *dimdef;
	const char *corename;

	char option_name[1024];
	char core_string[1024];
	char dim_args[2048];

	FILE *fd, *fd2, *fd3, *fd4, *fd5, *fd6;

	int in_subpool;

	corename = ezxml_attr(registry, "core");

	sprintf(core_string, "_%s_", corename);

	// For now, don't include core name in subroutines
	sprintf(core_string, "_");

	// Open files
	fd = fopen("read_dimensions.inc", "w+");
	fd2 = fopen("dim_dummy_args.inc", "w+");
	fd3 = fopen("add_dims_to_pool.inc", "w+");
	fd4 = fopen("dim_dummy_defines_input.inc", "w+");
	fd5 = fopen("dim_dummy_defines_noinput.inc", "w+");
	fd6 = fopen("dim_dummy_defines_inout.inc", "w+");

	dim_args[0] = '\0';

	// Parse dimensions that need to be read in
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			// Only dimensions that don't have a definition
			if(dimdef == NULL){
				// Write the read_dimension file, which reads dimensions from the input file
				fortprintf(fd, "      call mpas_io_inq_dim(inputHandle, '%s', %s, ierr)\n", dimname, dimname);

				// Write the dim dummy args list, which defines the argument list for subroutines that pass dimensions.
				if(strlen(dim_args) <= 1){
					sprintf(dim_args, "%s", dimname);
				} else {
					sprintf(dim_args, "%s, %s", dim_args, dimname);
				}


				fortprintf(fd3, "      call mpas_pool_add_dimension(block_ptr %% dimensions, '%s', %s)\n", dimname, dimname);

				// Write the dim dummy defines files, which defines the dimension within a subroutine.
				fortprintf(fd4, "      integer, intent(in) :: %s\n", dimname);
				fortprintf(fd5, "      integer :: %s\n", dimname);
				fortprintf(fd6, "      integer, intent(inout) :: %s\n", dimname);

			}
		}
	}

	fortprintf(fd2, "%s &\n", dim_args);

	// Close files
	fclose(fd);
	fclose(fd2);
	fclose(fd3);
	fclose(fd4);
	fclose(fd5);
	fclose(fd6);

	// Write subroutine to defined derived dimensions
	fd = fopen("define_dimensions.inc", "w+");
	fortprintf(fd, "   subroutine mpas_define%sderived_dimensions(dimensionPool, configPool)\n", core_string);
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: configPool\n");
	fortprintf(fd, "\n");

	// Define all dimensions
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			fortprintf(fd, "      integer, pointer :: %s\n", dimname);

			if(dimdef != NULL){
				// Namelist defined dimension
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);
					// Need to define a variable to hold the namelist value
					// First need to find the registry defined namlist option, so we can determine type:
					for (nmlrec_xml = ezxml_child(registry, "nml_record"); nmlrec_xml; nmlrec_xml = nmlrec_xml->next){
						nmlrecname = ezxml_attr(nmlrec_xml, "name");
						for (nmlopt_xml = ezxml_child(nmlrec_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
							nmloptname = ezxml_attr(nmlopt_xml, "name");
							nmlopttype = ezxml_attr(nmlopt_xml, "type");

							if(strncmp(option_name, nmloptname, 1024) == 0){
								if(strncmp(nmlopttype, "real", 1024) == 0){
									fortprintf(fd, "      real (kind=RKIND), pointer :: %s\n", nmloptname);
								} else if(strncmp(nmlopttype, "integer", 1024) == 0){
									fortprintf(fd, "      integer, pointer :: %s\n", nmloptname);
								} else if(strncmp(nmlopttype, "logical", 1024) == 0){
									fortprintf(fd, "      logical, pointer :: %s\n", nmloptname);
								} else if(strncmp(nmlopttype, "character", 1024) == 0){
									fortprintf(fd, "      character (len=StrKIND), pointer :: %s\n", nmloptname);
								}
							}
						}
					}
				}
			}
		}
	}

	fortprintf(fd,"\n");

	// Get values of all read in dimensions from pool
	// And config options from namelist defined dimensions.
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			if(dimdef == NULL){
				fortprintf(fd, "      call mpas_pool_get_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);
			} else {
				// Namelist defined dimension
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);
					// Need to define a variable to hold the namelist value
					// First need to find the registry defined namlist option, so we can determine type:
					for (nmlrec_xml = ezxml_child(registry, "nml_record"); nmlrec_xml; nmlrec_xml = nmlrec_xml->next){
						nmlrecname = ezxml_attr(nmlrec_xml, "name");
						nmlrecinsub = ezxml_attr(nmlrec_xml, "in_subpool");

						in_subpool = 0;

						if(nmlrecinsub && strcmp(nmlrecinsub, "true") == 0){
							in_subpool = 1;
						}

						for (nmlopt_xml = ezxml_child(nmlrec_xml, "nml_option"); nmlopt_xml; nmlopt_xml = nmlopt_xml->next){
							nmloptname = ezxml_attr(nmlopt_xml, "name");
							if(strcmp(option_name, nmloptname) == 0){
								if(in_subpool){
									fortprintf(fd, "      call mpas_pool_get_config(configPool, '%s', %s, '%s')\n", nmloptname, nmloptname, nmlrecname);
								} else {
									fortprintf(fd, "      call mpas_pool_get_config(configPool, '%s', %s)\n", nmloptname, nmloptname);
								}
							}
						}
					}
				}
			}
		}
	}

	fortprintf(fd,"\n");

	// Define and add dimensions to pool
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");
			dimdef = ezxml_attr(dim_xml, "definition");

			if(dimdef != NULL){
				fortprintf(fd, "      allocate(%s)\n", dimname);
				// Namelist defined dimension
				if(strncmp(dimdef, "namelist:", 9) == 0){
					snprintf(option_name, 1024, "%s", (dimdef)+9);
					fortprintf(fd, "      %s = %s\n", dimname, option_name);
				} else {
					fortprintf(fd, "      %s = %s\n", dimname, dimdef);
				}

				fortprintf(fd, "      call mpas_pool_add_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);
			}
		}
	}
	fortprintf(fd, "\n");
	fortprintf(fd, "   end subroutine mpas_define%sderived_dimensions\n", core_string);

	return 0;
}/*}}}*/


int parse_var_array(FILE *fd, ezxml_t registry, ezxml_t superStruct, ezxml_t varArray, const char * corename)/*{{{*/
{
	ezxml_t struct_xml, var_arr_xml, var_xml, var_xml2;
	ezxml_t packages_xml, package_xml;
	ezxml_t streams_xml, stream_xml, streams_xml2, stream_xml2;

	const char *structname, *structlevs, *structpackages;
	const char *substructname;
	const char *vararrname, *vararrtype, *vararrdims, *vararrpersistence, *vararrdefaultval, *vararrpackages;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *vardefaultval, *varpackages;
	const char *varname2, *vararrgroup2, *vararrname_in_code;
	const char *varname_in_code;
	const char *streamname, *streamname2;
	const char *packagename;
	const char *vararrtimelevs;

	int err;

	int iostreams;
	int time_lev, time_levs;
	int i, skip_var, skip_stream;
	int ndims, type, hasTime, decomp, in_stream;
	int persistence;
	char *string, *tofree, *token;
	char pointer_name[1024];
	char spacing[1024], sub_spacing[1024];
	char default_value[1024];

	structname = ezxml_attr(superStruct, "name");

	var_arr_xml = varArray;

	// All sub-structs have been parsed and generated at this point. Time to generate this struct
	// Start by generating variable arrays
	vararrname = ezxml_attr(var_arr_xml, "name");
	vararrtype = ezxml_attr(var_arr_xml, "type");
	vararrdims = ezxml_attr(var_arr_xml, "dimensions");
	vararrpersistence = ezxml_attr(var_arr_xml, "persistence");
	vararrdefaultval = ezxml_attr(var_arr_xml, "default_value");
	vararrpackages = ezxml_attr(var_arr_xml, "packages");
	vararrtimelevs = ezxml_attr(var_arr_xml, "time_levs");
	vararrname_in_code = ezxml_attr(var_arr_xml, "name_in_code");

	if(!vararrtimelevs){
		vararrtimelevs = ezxml_attr(superStruct, "time_levs");
	}

	if(vararrtimelevs){
		time_levs = atoi(vararrtimelevs);
		if(time_levs < 1){
			time_levs = 1;
		}
	} else {
		time_levs = 1;
	}

	if(!vararrname_in_code){
		vararrname_in_code = ezxml_attr(var_arr_xml, "name");
	}

	persistence = check_persistence(vararrpersistence);

	fortprintf(fd, "! Define var array %s\n", vararrname);
	snprintf(spacing, 1024, "      ");


	// Determine field type and default value.
	get_field_information(vararrtype, vararrdefaultval, default_value, &type);

	// Determine ndims, hasTime, and decomp type
	get_dimension_information(vararrdims, &ndims, &hasTime, &decomp);
	ndims++; // Add a dimension for constituents in var_array

	// Determine name of pointer for this field.
	set_pointer_name(type, ndims, pointer_name);
	fortprintf(fd, "      allocate(%s(%d))\n", pointer_name, time_levs);

	fortprintf(fd, "      index_counter = 0\n", spacing);
	fortprintf(fd, "      group_counter = -1\n", spacing);
	fortprintf(fd, "      group_start = -1\n", spacing);
	fortprintf(fd, "      group_started = .false.\n", spacing);
	fortprintf(fd, "\n");

	// Write index values and group counter values.
	// Define each array_group in contiguous sections.
	for (var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
		varname = ezxml_attr(var_xml, "name");
		varpackages = ezxml_attr(var_xml, "packages");
		vararrgroup = ezxml_attr(var_xml, "array_group");
		varname_in_code = ezxml_attr(var_xml, "name_in_code");
		skip_var = 0;

		if(!varname_in_code){
			varname_in_code = ezxml_attr(var_xml, "name");
		}

		for (var_xml2 = ezxml_child(var_arr_xml, "var"); var_xml2 && var_xml2 != var_xml; var_xml2 = var_xml2->next){
			// Check if the current array group has already been touched.
			vararrgroup2 = ezxml_attr(var_xml2, "array_group");

			if (strncmp(vararrgroup, vararrgroup2, 1024) == 0){
				skip_var = 1;
			}
		}

		if(!skip_var){

			fortprintf(fd, "! Starting group %s\n", vararrgroup);
			fortprintf(fd, "! Define constituent var %s\n", varname);
			fortprintf(fd, "! My Packages are %s\n", varpackages);

			// If no packages are defined, default to var_arr packages.
			if(varpackages == NULL){
				if(vararrpackages != NULL){
					varpackages = ezxml_attr(var_arr_xml, "packages");
				}
			}

			// Parse packages if they are defined
			sub_spacing[0] = '\0';
			if(varpackages){
				fortprintf(fd, "      if (");
				string = strdup(varpackages);
				tofree = string;
				token = strsep(&string, ";");
				fortprintf(fd, "%sActive", token);

				while( (token = strsep(&string, ";")) != NULL){
					fortprintf(fd, " .or. %sActive", token);
				}

				fortprintf(fd, ") then\n");
				snprintf(sub_spacing, 1024, "   ");
			}

			fortprintf(fd, "      %sindex_counter = index_counter + 1\n", sub_spacing);
			fortprintf(fd, "      %sif (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s   call mpas_pool_add_dimension(newSubPool, 'index_%s', index_counter)\n", sub_spacing, varname_in_code);
			fortprintf(fd, "      %send if\n", sub_spacing);
			fortprintf(fd, "      %sgroup_counter = group_counter + 1\n", sub_spacing);
			fortprintf(fd, "      %sif (.not. group_started) then\n", sub_spacing);
			fortprintf(fd, "      %s   group_start = index_counter\n", sub_spacing);
			fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_start', group_start)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s   end if\n", sub_spacing);
			fortprintf(fd, "      %s   group_started = .true.\n", sub_spacing);
			fortprintf(fd, "      %send if\n", sub_spacing);

			// If Packages are defined, write else clause
			if(varpackages){
				fortprintf(fd, "   %selse\n", sub_spacing);
				fortprintf(fd, "      %s  if (associated(newSubPool)) then\n", sub_spacing);
				fortprintf(fd, "      %s     call mpas_pool_add_dimension(newSubPool, 'index_%s', -1)\n", sub_spacing, varname_in_code);
				fortprintf(fd, "      %s  end if\n", sub_spacing);
				fortprintf(fd, "   %send if\n", sub_spacing);
			}

			// Add the rest of the variables from the current group.
			if(var_xml->next){
				for(var_xml2 = var_xml->next; var_xml2; var_xml2 = var_xml2->next){
					vararrgroup2 = ezxml_attr(var_xml2, "array_group");

					// var_xml2 is in the current array group
					if(strncmp(vararrgroup, vararrgroup2, 1024) == 0){
						varname = ezxml_attr(var_xml2, "name");
						varpackages = ezxml_attr(var_xml2, "packages");
						varname_in_code = ezxml_attr(var_xml2, "name_in_code");

						if(!varname_in_code){
							varname_in_code = ezxml_attr(var_xml2, "name");
						}


						// If no packages are defined, default to var_arr packages.
						if(varpackages == NULL){
							if(vararrpackages != NULL){
								varpackages = ezxml_attr(var_arr_xml, "packages");
							}
						}

						fortprintf(fd, "! Define constituent var %s\n", varname);
						fortprintf(fd, "! My packages are %s\n", varpackages);

						// Parse packages if they are defined
						sub_spacing[0] = '\0';
						if(varpackages){
							fortprintf(fd, "%sif (", spacing);
							string = strdup(varpackages);
							tofree = string;
							token = strsep(&string, ";");
							fortprintf(fd, "%sActive", token);

							while( (token = strsep(&string, ";")) != NULL){
								fortprintf(fd, " .or. %sActive", token);
							}

							fortprintf(fd, ") then\n");
							snprintf(sub_spacing, 1024, "   ");
						}

						fortprintf(fd, "      %sindex_counter = index_counter + 1\n", sub_spacing);
						fortprintf(fd, "      %sif (associated(newSubPool)) then\n", sub_spacing);
						fortprintf(fd, "      %s   call mpas_pool_add_dimension(newSubPool, 'index_%s', index_counter)\n", sub_spacing, varname_in_code);
						fortprintf(fd, "      %send if\n", sub_spacing);
						fortprintf(fd, "      %sgroup_counter = group_counter + 1\n", sub_spacing);
						fortprintf(fd, "      %sif (.not. group_started) then\n", sub_spacing);
						fortprintf(fd, "      %s   group_start = index_counter\n", sub_spacing);
						fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
						fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_start', group_start)\n", sub_spacing, vararrgroup);
						fortprintf(fd, "      %s   end if\n", sub_spacing);
						fortprintf(fd, "      %s   group_started = .true.\n", sub_spacing);
						fortprintf(fd, "      %send if\n", sub_spacing);

						// If Packages are defined, write else clause
						if(varpackages != NULL){
							fortprintf(fd, "   %selse\n", sub_spacing);
							fortprintf(fd, "   %s   if (associated(newSubPool)) then\n", sub_spacing);
							fortprintf(fd, "   %s      call mpas_pool_add_dimension(newSubPool, 'index_%s', -1)\n", sub_spacing, varname_in_code);
							fortprintf(fd, "   %s   end if\n", sub_spacing);
							fortprintf(fd, "   %send if\n", sub_spacing);
						}
					}
				}
			}

			fortprintf(fd, "      %sif (.not. group_started) then\n", sub_spacing);
			fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_start', -1)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_end', -1)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s   end if\n", sub_spacing);
			fortprintf(fd, "      %selse\n", sub_spacing);
			fortprintf(fd, "      %s   group_started = .false.\n", sub_spacing);
			fortprintf(fd, "      %s   if (associated(newSubPool)) then\n", sub_spacing);
			fortprintf(fd, "      %s      call mpas_pool_add_dimension(newSubPool, '%s_end', index_counter)\n", sub_spacing, vararrgroup);
			fortprintf(fd, "      %s   end if\n", sub_spacing);
			fortprintf(fd, "      %send if\n", sub_spacing);
			fortprintf(fd, "! End of group       \n", vararrgroup);
		}
	}

	fortprintf(fd, "\n");

	// Setup constituent names
	fortprintf(fd, "      numConstituents = index_counter\n");
	fortprintf(fd, "      if (associated(newSubPool)) then\n");
	fortprintf(fd, "         call mpas_pool_add_dimension(newSubPool, 'num_%s', numConstituents)\n", vararrname);
	fortprintf(fd, "      end if\n");

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		fortprintf(fd, "! Defining time level %d\n", time_lev);
		fortprintf(fd, "      allocate( %s(%d) %% constituentNames(numConstituents) )\n", pointer_name, time_lev);
		fortprintf(fd, "      allocate(%s(%d) %% ioinfo)\n", pointer_name, time_lev);
		fortprintf(fd, "      %s(%d) %% fieldName = '%s'\n", pointer_name, time_lev, vararrname);
		if (hasTime) {
			fortprintf(fd, "      %s(%d) %% hasTimeDimension = .true.\n", pointer_name, time_lev);
		} else {
			fortprintf(fd, "      %s(%d) %% hasTimeDimension = .false.\n", pointer_name, time_lev);
		}
		fortprintf(fd, "      %s(%d) %% isVarArray = .true.\n", pointer_name, time_lev);
		if(ndims > 0){
			if(persistence == SCRATCH){
				fortprintf(fd, "      %s(%d) %% isPersistent = .false.\n", pointer_name, time_lev);
				fortprintf(fd, "      %s(%d) %% isActive = .false.\n", pointer_name, time_lev);
			} else {
				fortprintf(fd, "      %s(%d) %% isPersistent = .true.\n", pointer_name, time_lev);
				fortprintf(fd, "      %s(%d) %% isActive = .false.\n", pointer_name, time_lev);
			}
		}
		fortprintf(fd, "\n");
		for(var_xml = ezxml_child(var_arr_xml, "var"); var_xml; var_xml = var_xml->next){
			varname = ezxml_attr(var_xml, "name");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");

			if(!varname_in_code){
				varname_in_code = ezxml_attr(var_xml, "name");
			}

			fortprintf(fd, "      if (associated(newSubPool)) then\n");
			fortprintf(fd, "         call mpas_pool_get_dimension(newSubPool, 'index_%s', const_index)\n", varname_in_code);
			fortprintf(fd, "      end if\n");
			fortprintf(fd, "      if (index_counter > 0) then\n", spacing);
			fortprintf(fd, "         %s(%d) %% constituentNames(const_index) = '%s'\n", pointer_name, time_lev, varname);
			fortprintf(fd, "      end if\n", spacing);
		}

		fortprintf(fd, "\n");

		// Setup dimensions
		fortprintf(fd, "! Setup dimensions for       \n", vararrname);
		i = 1;
		fortprintf(fd, "      %s(%d) %% dimNames(%d) = 'num_%s'\n", pointer_name, time_lev, i, vararrname);
		fortprintf(fd, "      %s(%d) %% dimSizes(%d) = numConstituents\n", pointer_name, time_lev, i);

		string = strdup(vararrdims);
		tofree = string;
		token = strsep(&string, " ");

		if(strncmp(token, "Time", 1024) != 0){
			i++;
			if(strncmp(token, "nCells", 1024) == 0 || strncmp(token, "nEdges", 1024) == 0 || strncmp(token, "nVertices", 1024) == 0){
				fortprintf(fd, "      %s(%d) %% dimNames(%d) = '%s'\n", pointer_name, time_lev, i, token);
				fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s+1\n", pointer_name, time_lev, i, token);
			} else {
				fortprintf(fd, "      %s(%d) %% dimNames(%d) = '%s'\n", pointer_name, time_lev, i, token);
				fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s\n", pointer_name, time_lev, i, token);
			}
		}
		while( (token = strsep(&string, " ")) != NULL){
			if(strncmp(token, "Time", 1024) != 0){
				i++;
				if(strncmp(token, "nCells", 1024) == 0 || strncmp(token, "nEdges", 1024) == 0 || strncmp(token, "nVertices", 1024) == 0){
					fortprintf(fd, "      %s(%d) %% dimNames(%d) = '%s'\n", pointer_name, time_lev, i, token);
					fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s+1\n", pointer_name, time_lev, i, token);
				} else {
					fortprintf(fd, "      %s(%d) %% dimNames(%d) = '%s'\n", pointer_name, time_lev, i, token);
					fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s\n", pointer_name, time_lev, i, token);
				}
			}
		}
		free(tofree);

		fortprintf(fd, "\n");

		// Setup array pointer
		fortprintf(fd, "! Allocate space for data\n");
		if(!vararrpersistence || strcmp(vararrpersistence, "scratch") != 0){
			switch(ndims){
				default:
					break;
				case 1:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1)))\n", pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 2:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2)))\n", pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 3:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2), %s(%d) %% dimSizes(3)))\n", pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 4:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2), %s(%d) %% dimSizes(3), %s(%d) %% dimSizes(4)))\n", pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 5:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2), %s(%d) %% dimSizes(3), %s(%d) %% dimSizes(4), %s(%d) %% dimSizes(5)))\n", pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
			}

			fortprintf(fd, "      %s(%d) %% array = %s\n", pointer_name, time_lev, default_value);

		} else {
			if(ndims > 0){
				fortprintf(fd, "      nullify(%s(%d) %% array)\n", pointer_name, time_lev);
			}
		}

		fortprintf(fd, "      nullify(%s(%d) %% next)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% prev)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% sendList)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% recvList)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% copyList)\n", pointer_name, time_lev);
		fortprintf(fd, "      %s(%d) %% block => block\n", pointer_name, time_lev);
	}

	// Parse packages if they are defined
	fortprintf(fd, "\n");
	spacing[0] = '\0';
	if(vararrpackages != NULL){
		fortprintf(fd, "      if (");
		string = strdup(vararrpackages);
		tofree = string;
		token = strsep(&string, ";");
		fortprintf(fd, "%sActive", token);

		while( (token = strsep(&string, ";")) != NULL){
			fortprintf(fd, " .or. %sActive", token);
		}

		fortprintf(fd, ") then\n");
		snprintf(spacing, 1024, "         ");
	}


	// Add field to pool
	fortprintf(fd, "! Add field to pool\n");
	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		fortprintf(fd, "%s%s(%d) %% isActive = .true.\n", spacing, pointer_name, time_lev);
	}
	fortprintf(fd, "%scall mpas_pool_add_field(newSubPool, '%s', %s)\n", spacing, vararrname_in_code, pointer_name);

	if(vararrpackages != NULL) {
		fortprintf(fd, "      end if\n");
	}

	fortprintf(fd, "      call mpas_pool_add_field(block %% allFields, '%s', %s)\n", vararrname, pointer_name);
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int parse_var(FILE *fd, ezxml_t registry, ezxml_t superStruct, ezxml_t currentVar, const char * corename)/*{{{*/
{
	ezxml_t struct_xml, var_xml, var_xml2;
	ezxml_t packages_xml, package_xml;
	ezxml_t streams_xml, stream_xml, streams_xml2, stream_xml2;

	const char *structtimelevs, *vartimelevs;
	const char *structname, *structlevs, *structpackages;
	const char *substructname;
	const char *varname, *varpersistence, *vartype, *vardims, *varunits, *vardesc, *vararrgroup, *varstreams, *vardefaultval, *varpackages;
	const char *varname2, *vararrgroup2;
	const char *varname_in_code;
	const char *streamname, *streamname2;
	const char *packagename;

	int err;

	int iostreams;
	int i, skip_var, skip_stream;
	int time_lev, time_levs;
	int ndims, type, hasTime, decomp, in_stream;
	int persistence;
	char *string, *tofree, *token;
	char pointer_name[1024];
	char package_spacing[1024];
	char default_value[1024];

	var_xml = currentVar;

	structname = ezxml_attr(superStruct, "name");
	structtimelevs = ezxml_attr(superStruct, "time_levs");

	// Define independent variables
	varname = ezxml_attr(var_xml, "name");
	vartype = ezxml_attr(var_xml, "type");
	vardims = ezxml_attr(var_xml, "dimensions");
	varpersistence = ezxml_attr(var_xml, "persistence");
	varpackages = ezxml_attr(var_xml, "packages");
	vardefaultval = ezxml_attr(var_xml, "default_value");
	vartimelevs = ezxml_attr(var_xml, "time_levs");
	varname_in_code = ezxml_attr(var_xml, "name_in_code");

	if(!varname_in_code){
		varname_in_code = ezxml_attr(var_xml, "name");
	}

	if(!vartimelevs){
		vartimelevs = ezxml_attr(superStruct, "time_levs");
	}

	if(vartimelevs){
		time_levs = atoi(vartimelevs);
		if(time_levs < 1){
			time_levs = 1;
		}
	} else {
		time_levs = 1;
	}

	persistence = check_persistence(varpersistence);

	fortprintf(fd, "! Define variable %s\n", varname);


	// Determine field type and default value.
	get_field_information(vartype, vardefaultval, default_value, &type);

	// Determine ndims, hasTime, and decomp type
	get_dimension_information(vardims, &ndims, &hasTime, &decomp);

	// Determine name of pointer for this field.
	set_pointer_name(type, ndims, pointer_name);
	fortprintf(fd, "      allocate(%s(%d))\n", pointer_name, time_levs);

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		fortprintf(fd, "\n");
		fortprintf(fd, "! Setting up time level %d\n", time_lev);
		fortprintf(fd, "      allocate(%s(%d) %% ioinfo)\n",  pointer_name, time_lev);
		fortprintf(fd, "      %s(%d) %% fieldName = '%s'\n", pointer_name, time_lev, varname);
		fortprintf(fd, "      %s(%d) %% isVarArray = .false.\n", pointer_name, time_lev);
		if(hasTime) {
			fortprintf(fd, "      %s(%d) %% hasTimeDimension = .true.\n", pointer_name, time_lev);
		} else {
			fortprintf(fd, "      %s(%d) %% hasTimeDimension = .false.\n", pointer_name, time_lev);
		}

		if(ndims > 0){
			if(persistence == SCRATCH){
				fortprintf(fd, "      %s(%d) %% isPersistent = .false.\n", pointer_name, time_lev);
				fortprintf(fd, "      %s(%d) %% isActive = .false.\n", pointer_name, time_lev);
			} else {
				fortprintf(fd, "      %s(%d) %% isPersistent = .true.\n", pointer_name, time_lev);
				fortprintf(fd, "      %s(%d) %% isActive = .false.\n", pointer_name, time_lev);
			}

			// Setup dimensions
			fortprintf(fd, "! Setting up dimensions\n");
			string = strdup(vardims);
			tofree = string;
			i = 1;
			token = strsep(&string, " ");
			if(strncmp(token, "Time", 1024) != 0){
				fortprintf(fd, "      %s(%d) %% dimNames(%d) = '%s'\n", pointer_name, time_lev, i, token);
				if(strcmp(token, "nCells") == 0 || strcmp(token, "nEdges") == 0 || strcmp(token, "nVertices") == 0){
					fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s + 1\n", pointer_name, time_lev, i, token);
				} else {
					fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s\n", pointer_name, time_lev, i, token);
				}
				i++;
			}
			while( (token = strsep(&string, " ")) != NULL){
				if(strncmp(token, "Time", 1024) != 0){
					fortprintf(fd, "      %s(%d) %% dimNames(%d) = '%s'\n", pointer_name, time_lev, i, token);
					if(strcmp(token, "nCells") == 0 || strcmp(token, "nEdges") == 0 || strcmp(token, "nVertices") == 0){
						fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s + 1\n", pointer_name, time_lev, i, token);
					} else {
						fortprintf(fd, "      %s(%d) %% dimSizes(%d) = %s\n", pointer_name, time_lev, i, token);
					}
					i++;
				}
			}
			free(tofree);



			fortprintf(fd, "! Allocate space for data\n");
			switch(ndims){
				default:
					break;
				case 1:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1)))\n", pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 2:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2)))\n",
							pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 3:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2), %s(%d) %% dimSizes(3)))\n",
							pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 4:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2), %s(%d) %% dimSizes(3), %s(%d) %% dimSizes(4)))\n",
							pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
				case 5:
					fortprintf(fd, "      allocate(%s(%d) %% array(%s(%d) %% dimSizes(1), %s(%d) %% dimSizes(2), %s(%d) %% dimSizes(3), %s(%d) %% dimSizes(4), %s(%d) %% dimSizes(5)))\n",
							pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev, pointer_name, time_lev);
					break;
			}

			fortprintf(fd, "      %s(%d) %% array = %s\n", pointer_name, time_lev, default_value);
		} else if(ndims == 0){
			fortprintf(fd, "      %s(%d) %% scalar = %s\n", pointer_name, time_lev, default_value);
		}


		fortprintf(fd, "      nullify(%s(%d) %% next)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% prev)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% sendList)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% recvList)\n", pointer_name, time_lev);
		fortprintf(fd, "      nullify(%s(%d) %% copyList)\n", pointer_name, time_lev);
		fortprintf(fd, "      %s(%d) %% block => block\n", pointer_name, time_lev);
	}

	// Parse packages if they are defined
	fortprintf(fd, "\n");
	snprintf(package_spacing, 1024, "      ");
	if(varpackages != NULL){
		fortprintf(fd, "      if (");
		string = strdup(varpackages);
		tofree = string;
		token = strsep(&string, ";");
		fortprintf(fd, "%sActive", token);

		while( (token = strsep(&string, ";")) != NULL){
			fortprintf(fd, " .or. %sActive", token);
		}

		fortprintf(fd, ") then\n");
		snprintf(package_spacing, 1024, "         ");
	}

	for(time_lev = 1; time_lev <= time_levs; time_lev++){
		fortprintf(fd, "%s%s(%d) %% isActive = .true.\n", package_spacing, pointer_name, time_lev);
	}
	fortprintf(fd, "%scall mpas_pool_add_field(newSubPool, '%s', %s)\n", package_spacing , varname_in_code, pointer_name);

	if(varpackages != NULL){
		fortprintf(fd, "      end if\n");
	}
	fortprintf(fd, "      call mpas_pool_add_field(block %% allFields, '%s', %s)\n", varname, pointer_name);
	fortprintf(fd, "\n");

	return 0;
}/*}}}*/


int parse_struct(FILE *fd, ezxml_t registry, ezxml_t superStruct, int subpool, const char *parentname, const char * corename)/*{{{*/
{
	ezxml_t dims_xml, dim_xml;
	ezxml_t struct_xml, var_arr_xml, var_xml, var_xml2;
	ezxml_t struct_xml2;
	ezxml_t packages_xml, package_xml;

	const char *dimname;
	const char *structname, *structlevs, *structpackages;
	const char *structname2;
	const char *substructname;
	const char *streamname, *streamname2;
	const char *packagename;
	const char *structnameincode;

	char *string, *tofree, *token;
	char spacing[1024];
	char core_string[1024];
	char pool_name[1024];
	char package_list[2048];

	int skip_struct, no_packages;
	int err;

	sprintf(core_string, "_%s_", corename);

	// For now, don't include core name in subroutines
	sprintf(core_string, "_");

	if(subpool){
		sprintf(pool_name, "%s_subpool", parentname);
	} else {
		sprintf(pool_name, "pool");
	}

	structname = ezxml_attr(superStruct, "name");
	structnameincode = ezxml_attr(superStruct, "name_in_code");
	
	if(!structnameincode){
		structnameincode = ezxml_attr(superStruct, "name");
	}

	structpackages = ezxml_attr(superStruct, "packages");

	// Extract all sub structs
	for (struct_xml = ezxml_child(superStruct, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		err = parse_struct(fd, registry, struct_xml, 1, structname, corename);
	}

	fortprintf(fd, "   subroutine mpas_generate%s%s_%s(block, structPool, dimensionPool, packagePool)\n", core_string, pool_name, structname);
	fortprintf(fd, "      type (block_type), intent(inout), pointer :: block\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: structPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(in) :: packagePool\n");
	write_field_pointer_arrays(fd);
	fortprintf(fd, "      type (mpas_pool_type), pointer :: newSubPool\n");
	fortprintf(fd, "      type (mpas_pool_iterator_type) :: dimItr\n");
	fortprintf(fd, "      integer, pointer :: dim0D\n");
	fortprintf(fd, "      integer, dimension(:), pointer :: dim1D\n");
	fortprintf(fd, "      integer :: group_counter\n");
	fortprintf(fd, "      logical :: group_started\n");
	fortprintf(fd, "      integer :: group_start\n");
	fortprintf(fd, "      integer :: index_counter\n");
	fortprintf(fd, "      integer, pointer :: const_index\n");
	fortprintf(fd, "\n");

	// Need to define logicals for all packages
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");

			fortprintf(fd, "      logical, pointer :: %sActive\n", packagename);
		}
	}

	fortprintf(fd, "\n");

	// Need to define integers for all dimensions
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");

			fortprintf(fd, "      integer, pointer :: %s\n", dimname);
		}
	}

	fortprintf(fd, "\n");
	fortprintf(fd, "      integer :: numConstituents\n");
	fortprintf(fd, "\n");

	fortprintf(fd, "      nullify(newSubPool)\n");

	fortprintf(fd, "      group_counter = -1\n");
	fortprintf(fd, "      group_started = .false.\n");
	fortprintf(fd, "      group_start = -1\n");

	// Need to get value of package flags
	for (packages_xml = ezxml_child(registry, "packages"); packages_xml; packages_xml = packages_xml->next){
		for (package_xml = ezxml_child(packages_xml, "package"); package_xml; package_xml = package_xml->next){
			packagename = ezxml_attr(package_xml, "name");

			fortprintf(fd, "      call mpas_pool_get_package(packagePool, '%sActive', %sActive)\n", packagename, packagename);
		}
	}

	fortprintf(fd, "\n");

	// Parse packages if they are defined
	package_list[0] = '\0';
	no_packages = build_struct_package_lists(superStruct, package_list);

	spacing[0] = '\0';
	if(!no_packages){
		fortprintf(fd, "      if (");
		string = strdup(package_list);
		tofree = string;
		token = strsep(&string, ";");
		fortprintf(fd, "%sActive", token);

		while( (token = strsep(&string, ";")) != NULL){
			fortprintf(fd, " .or. %sActive", token);
		}

		fortprintf(fd, ") then\n");
		sprintf(spacing, "   ");
	}

	// Setup new pool to be added into structPool
	fortprintf(fd, "      %sallocate(newSubPool)\n", spacing);
	fortprintf(fd, "      %scall mpas_pool_create_pool(newSubPool)\n", spacing);
	fortprintf(fd, "      %scall mpas_pool_add_subpool(structPool, '%s', newSubPool)\n", spacing, structnameincode);
	fortprintf(fd, "      %scall mpas_pool_add_subpool(block %% allStructs, '%s', newSubPool)\n", spacing, structname);

	if(!no_packages){
		fortprintf(fd, "      end if\n");
	}

	// Need to get value of dimensions
	for (dims_xml = ezxml_child(registry, "dims"); dims_xml; dims_xml = dims_xml->next){
		for (dim_xml = ezxml_child(dims_xml, "dim"); dim_xml; dim_xml = dim_xml->next){
			dimname = ezxml_attr(dim_xml, "name");

			fortprintf(fd, "      call mpas_pool_get_dimension(dimensionPool, '%s', %s)\n", dimname, dimname);
		}
	}

	fortprintf(fd, "\n");

	// All sub-structs have been parsed and generated at this point. Time to generate this struct
	// Start by generating variable arrays
	for (var_arr_xml = ezxml_child(superStruct, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){
		parse_var_array(fd, registry, superStruct, var_arr_xml, corename);
	}


	// Define independent variables
	for (var_xml = ezxml_child(superStruct, "var"); var_xml; var_xml = var_xml->next){
		parse_var(fd, registry, superStruct, var_xml, corename);
	}

	fortprintf(fd, "\n");

	// Extract all sub structs
	for (struct_xml = ezxml_child(superStruct, "var_struct"); struct_xml; struct_xml = struct_xml->next){
		fortprintf(fd, "      call mpas_generate%s%s_subpool_%s(block, newSubPool, dimensionPool, packagePool)\n", core_string, structname, substructname);
	}

	fortprintf(fd, "\n");
	fortprintf(fd, "      if (associated(newSubPool)) then\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block %% domain %% on_a_sphere)\n");
	fortprintf(fd, "         call mpas_pool_add_config(newSubPool, 'sphere_radius', block %% domain %% sphere_radius)\n");
	fortprintf(fd, "         call mpas_pool_begin_iteration(dimensionPool)\n");
	fortprintf(fd, "         do while( mpas_pool_get_next_member(dimensionPool, dimItr) )\n");
	fortprintf(fd, "            if (dimItr %% memberType == MPAS_POOL_DIMENSION) then\n");
	fortprintf(fd, "               if (dimItr %% nDims == 0) then\n");
	fortprintf(fd, "                  call mpas_pool_get_dimension(dimensionPool, dimItr %% memberName, dim0d)\n");
	fortprintf(fd, "                  call mpas_pool_add_dimension(newSubPool, dimItr %% memberName, dim0d)\n");
	fortprintf(fd, "               else if (dimItr %% nDims == 1) then\n");
	fortprintf(fd, "                  call mpas_pool_get_dimension(dimensionPool, dimItr %% memberName, dim1d)\n");
	fortprintf(fd, "                  call mpas_pool_add_dimension(newSubPool, dimItr %% memberName, dim1d)\n");
	fortprintf(fd, "               end if\n");
	fortprintf(fd, "            end if\n");
	fortprintf(fd, "         end do\n");
	fortprintf(fd, "      end if\n");
	fortprintf(fd, "\n");

	fortprintf(fd, "   end subroutine mpas_generate%s%s_%s\n", core_string, pool_name, structname);
	fortprintf(fd, "\n\n");

	return 0;
}/*}}}*/


int determine_struct_depth(int curLevel, ezxml_t superStruct){/*{{{*/
	ezxml_t subStruct;
	int max_depth, depth;

	max_depth = curLevel;

	for(subStruct = ezxml_child(superStruct, "var_struct"); subStruct; subStruct = subStruct->next){
		depth = determine_struct_depth(curLevel+1, subStruct);

		if(depth > max_depth){
			max_depth = depth;
		}
	}

	return max_depth;
}/*}}}*/


int generate_struct_links(FILE *fd, int curLevel, ezxml_t superStruct){/*{{{*/
	ezxml_t subStruct;
	ezxml_t var_arr_xml, var_xml;
	const char *structname;
	const char *vartimelevs;
	const char *varname, *vardims, *vartype;
	const char *vardefaultval;
	const char *varname_in_code;
	int depth;
	int err;
	int has_time;
	int time_lev, time_levs;
	int ndims, type;
	int decomp;
	char *string, *tofree, *token;
	char pointer_name[1024];
	char default_value[1024];

	depth = curLevel + 1;

	for(subStruct = ezxml_child(superStruct, "var_struct"); subStruct; subStruct = subStruct->next){
		structname = ezxml_attr(subStruct, "name");
		fortprintf(fd, "! ----------- NEW STRUCT ---------\n");
		fortprintf(fd, "! Get pointers to pools for struct %s\n", structname);
		fortprintf(fd, "! --------------------------------\n");
		if(curLevel == 0){
			fortprintf(fd, "      call mpas_pool_get_subpool(currentBlock %% structs, '%s', poolLevel%d)\n", structname, curLevel+1);
			fortprintf(fd, "      if(associated(prevBlock)) then\n");
			fortprintf(fd, "         call mpas_pool_get_subpool(prevBlock %% structs, '%s', prevPoolLevel%d)\n", structname, curLevel+1);
			fortprintf(fd, "      else\n");
			fortprintf(fd, "         nullify(prevPoolLevel%d)\n", curLevel+1);
			fortprintf(fd, "      end if\n");
			fortprintf(fd, "      if(associated(nextBlock)) then\n");
			fortprintf(fd, "         call mpas_pool_get_subpool(nextBlock %% structs, '%s', nextPoolLevel%d)\n", structname, curLevel+1);
			fortprintf(fd, "      else\n");
			fortprintf(fd, "         nullify(nextPoolLevel%d)\n", curLevel+1);
			fortprintf(fd, "      end if\n");
		} else {
			fortprintf(fd, "      call mpas_pool_get_subpool(poolLevel%d, '%s', poolLevel%d)\n", curLevel, structname, curLevel+1);
			fortprintf(fd, "      if(associated(prevBlock)) then\n");
			fortprintf(fd, "         call mpas_pool_get_subpool(prevPoolLevel%d, '%s', prevPoolLevel%d)\n", curLevel, structname, curLevel+1);
			fortprintf(fd, "      else\n");
			fortprintf(fd, "         nullify(prevPoolLevel%d)\n", curLevel+1);
			fortprintf(fd, "      end if\n");
			fortprintf(fd, "      if(associated(nextBlock)) then\n");
			fortprintf(fd, "         call mpas_pool_get_subpool(nextPoolLevel%d, '%s', nextPoolLevel%d)\n", curLevel, structname, curLevel+1);
			fortprintf(fd, "      else\n");
			fortprintf(fd, "         nullify(nextPoolLevel%d)\n", curLevel+1);
			fortprintf(fd, "      end if\n");
		}

		fortprintf(fd, "\n");
		// Link var arrays
		for(var_arr_xml = ezxml_child(subStruct, "var_array"); var_arr_xml; var_arr_xml = var_arr_xml->next){/*{{{*/
			varname = ezxml_attr(var_arr_xml, "name");
			vardims = ezxml_attr(var_arr_xml, "dimensions");
			vartimelevs = ezxml_attr(var_arr_xml, "time_levs");
			vartype = ezxml_attr(var_arr_xml, "type");
			vardefaultval = ezxml_attr(var_arr_xml, "default_value");

			if(!vartimelevs){
				vartimelevs = ezxml_attr(subStruct, "time_levs");
			}

			if(vartimelevs){
				time_levs = atoi(vartimelevs);
				if(time_levs < 1){
					time_levs = 1;
				}
			} else {
				time_levs = 1;
			}

			// Determine field type and default value.
			get_field_information(vartype, vardefaultval, default_value, &type);

			// Determine number of dimensions
			// and decomp type
			get_dimension_information(vardims, &ndims, &has_time, &decomp);
			ndims++; // Add a dimension for var_arrays

			// Using type and ndims, determine name of pointer for field.
			set_pointer_name(type, ndims, pointer_name);

			for(time_lev = 1; time_lev <= time_levs; time_lev++){
				fortprintf(fd, "! Linking %s for time level %d\n", varname, time_lev);
				fortprintf(fd, "      call mpas_pool_get_field(poolLevel%d, '%s', %s, %d)\n", curLevel+1, varname, pointer_name, time_lev);
				fortprintf(fd, "      if(associated(%s)) then\n", pointer_name);
				fortprintf(fd, "#ifdef MPAS_DEBUG\n");
				fortprintf(fd, "         write(stderrUnit,*) 'Linking %s'\n", varname);
				fortprintf(fd, "#endif\n");
				fortprintf(fd, "         if(associated(prevBlock)) then\n");
				fortprintf(fd, "            call mpas_pool_get_field(prevPoolLevel%d, '%s', %s %% prev, %d)\n", curLevel+1, varname, pointer_name, time_lev);
				fortprintf(fd, "         end if\n");
				fortprintf(fd, "         if(associated(nextBlock)) then\n");
				fortprintf(fd, "            call mpas_pool_get_field(nextPoolLevel%d, '%s', %s %% next, %d)\n", curLevel+1, varname, pointer_name, time_lev);
				fortprintf(fd, "         end if\n");

				if(decomp == CELLS){
					fortprintf(fd, "         %s %% sendList => currentBlock %% parinfo %% cellsToSend\n", pointer_name);
					fortprintf(fd, "         %s %% recvList => currentBlock %% parinfo %% cellsToRecv\n", pointer_name);
					fortprintf(fd, "         %s %% copyList => currentBlock %% parinfo %% cellsToCopy\n", pointer_name);
				} else if(decomp == EDGES){
					fortprintf(fd, "         %s %% sendList => currentBlock %% parinfo %% edgesToSend\n", pointer_name);
					fortprintf(fd, "         %s %% recvList => currentBlock %% parinfo %% edgesToRecv\n", pointer_name);
					fortprintf(fd, "         %s %% copyList => currentBlock %% parinfo %% edgesToCopy\n", pointer_name);
				} else if(decomp == VERTICES){
					fortprintf(fd, "         %s %% sendList => currentBlock %% parinfo %% verticesToSend\n", pointer_name);
					fortprintf(fd, "         %s %% recvList => currentBlock %% parinfo %% verticesToRecv\n", pointer_name);
					fortprintf(fd, "         %s %% copyList => currentBlock %% parinfo %% verticesToCopy\n", pointer_name);
				}

				fortprintf(fd, "      end if\n");
			}

			fortprintf(fd, "\n");
		}/*}}}*/

		// Link independent vars
		for(var_xml = ezxml_child(subStruct, "var"); var_xml; var_xml = var_xml->next){/*{{{*/
			varname = ezxml_attr(var_xml, "name");
			vardims = ezxml_attr(var_xml, "dimensions");
			vartimelevs = ezxml_attr(var_xml, "time_levs");
			vartype = ezxml_attr(var_xml, "type");
			vardefaultval = ezxml_attr(var_xml, "default_value");
			varname_in_code = ezxml_attr(var_xml, "name_in_code");

			if(!vartimelevs){
				vartimelevs = ezxml_attr(subStruct, "time_levs");
			}

			if(vartimelevs){
				time_levs = atoi(vartimelevs);
				if(time_levs < 1){
					time_levs = 1;
				}
			} else {
				time_levs = 1;
			}

			if(!varname_in_code){
				varname_in_code = ezxml_attr(var_xml, "name");
			}

			// Determine field type and default value.
			get_field_information(vartype, vardefaultval, default_value, &type);

			// Determine number of dimensions
			// and decomp type
			get_dimension_information(vardims, &ndims, &has_time, &decomp);

			// Using type and ndims, determine name of pointer for field.
			set_pointer_name(type, ndims, pointer_name);

			for(time_lev = 1; time_lev <= time_levs; time_lev++){
				fortprintf(fd, "! Linking %s for time level %d with name\n", varname, time_lev, varname_in_code);
				fortprintf(fd, "#ifdef MPAS_DEBUG\n");
				fortprintf(fd, "      write(stderrUnit,*) 'Linking %s with name %s'\n", varname, varname_in_code);
				fortprintf(fd, "#endif\n");
				fortprintf(fd, "      call mpas_pool_get_field(poolLevel%d, '%s', %s, %d)\n", curLevel+1, varname_in_code, pointer_name, time_lev);
				fortprintf(fd, "      if(associated(%s)) then\n", pointer_name);
				fortprintf(fd, "         if(associated(prevBlock)) then\n");
				fortprintf(fd, "            call mpas_pool_get_field(prevPoolLevel%d, '%s', %s %% prev, %d)\n", curLevel+1, varname_in_code, pointer_name, time_lev);
				fortprintf(fd, "         end if\n");
				fortprintf(fd, "         if(associated(nextBlock)) then\n");
				fortprintf(fd, "            call mpas_pool_get_field(nextPoolLevel%d, '%s', %s %% next, %d)\n", curLevel+1, varname_in_code, pointer_name, time_lev);
				fortprintf(fd, "         end if\n");

				if(decomp == CELLS){
					fortprintf(fd, "         %s %% sendList => currentBlock %% parinfo %% cellsToSend\n", pointer_name);
					fortprintf(fd, "         %s %% recvList => currentBlock %% parinfo %% cellsToRecv\n", pointer_name);
					fortprintf(fd, "         %s %% copyList => currentBlock %% parinfo %% cellsToCopy\n", pointer_name);
				} else if(decomp == EDGES){
					fortprintf(fd, "         %s %% sendList => currentBlock %% parinfo %% edgesToSend\n", pointer_name);
					fortprintf(fd, "         %s %% recvList => currentBlock %% parinfo %% edgesToRecv\n", pointer_name);
					fortprintf(fd, "         %s %% copyList => currentBlock %% parinfo %% edgesToCopy\n", pointer_name);
				} else if(decomp == VERTICES){
					fortprintf(fd, "         %s %% sendList => currentBlock %% parinfo %% verticesToSend\n", pointer_name);
					fortprintf(fd, "         %s %% recvList => currentBlock %% parinfo %% verticesToRecv\n", pointer_name);
					fortprintf(fd, "         %s %% copyList => currentBlock %% parinfo %% verticesToCopy\n", pointer_name);
				}
				fortprintf(fd, "      end if\n");

				fortprintf(fd, "\n");
			}
		}/*}}}*/

		err = generate_struct_links(fd, curLevel+1, subStruct);
	}

	return 0;
}/*}}}*/


int generate_field_links(ezxml_t registry){/*{{{*/
	ezxml_t struct_xml;
	const char *corename;
	FILE *fd;
	int i, structDepth, err;

	char core_string[1024];

	structDepth = determine_struct_depth(0, registry);

	corename = ezxml_attr(registry, "core");

	sprintf(core_string, "_%s_", corename);

	// For now, don't include core name in subroutines
	sprintf(core_string, "_");

	fd = fopen("link_fields.inc", "w+");

	fortprintf(fd, "   subroutine mpas%slink_fields(domain)\n", core_string);
	fortprintf(fd, "      type (domain_type), intent(in) :: domain\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      type (block_type), pointer :: blockPtr\n");
	fortprintf(fd, "      type (block_type), pointer :: prevBlock\n");
	fortprintf(fd, "      type (block_type), pointer :: nextBlock\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "      blockPtr => domain %% blocklist\n");
	fortprintf(fd, "      do while(associated(blockPtr))\n");
	fortprintf(fd, "         if (associated(blockPtr %% prev)) then\n");
	fortprintf(fd, "            prevBlock => blockPtr %% prev\n");
	fortprintf(fd, "         else\n");
	fortprintf(fd, "            nullify(prevBlock)\n");
	fortprintf(fd, "         end if\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "         if (associated(blockPtr %% next)) then\n");
	fortprintf(fd, "            nextBlock => blockPtr %% next\n");
	fortprintf(fd, "         else\n");
	fortprintf(fd, "            nullify(nextBlock)\n");
	fortprintf(fd, "         end if\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "         call mpas%slink_blocks(blockPtr, prevBlock, nextBlock)\n", core_string);
	fortprintf(fd, "\n");
	fortprintf(fd, "         blockPtr => blockPtr %% next\n");
	fortprintf(fd, "      end do\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "   end subroutine mpas%slink_fields\n", core_string);
	fortprintf(fd, "\n");
	fortprintf(fd, "\n");
	fortprintf(fd, "   subroutine mpas%slink_blocks(currentBlock, prevBlock, nextBlock)\n", core_string);
	fortprintf(fd, "      type (block_type), pointer, intent(in) :: currentBlock\n");
	fortprintf(fd, "      type (block_type), pointer, intent(in) :: prevBlock\n");
	fortprintf(fd, "      type (block_type), pointer, intent(in) :: nextBlock\n");
	write_field_pointers(fd);
	for(i = 1; i <= structDepth; i++){
		fortprintf(fd, "      type (mpas_pool_type), pointer :: poolLevel%d\n", i);
		fortprintf(fd, "      type (mpas_pool_type), pointer :: prevPoolLevel%d\n", i);
		fortprintf(fd, "      type (mpas_pool_type), pointer :: nextPoolLevel%d\n", i);
	}
	fortprintf(fd, "\n");

	err = generate_struct_links(fd, 0, registry);

	fortprintf(fd, "   end subroutine mpas%slink_blocks\n", core_string);

	fclose(fd);

	return 0;
}/*}}}*/


int generate_immutable_struct_contents(FILE *fd, const char *streamname, ezxml_t varstruct_xml){/*{{{*/
	ezxml_t var_xml, vararr_xml, substruct_xml;

	const char *optname, *optstream;

	/* Loop over fields looking for any that belong to the stream */
	for (vararr_xml = ezxml_child(varstruct_xml, "var"); vararr_xml; vararr_xml = vararr_xml->next) {
		optstream = ezxml_attr(vararr_xml, "streams");
		if (optstream != NULL && strstr(optstream, streamname) != NULL) {
			optname = ezxml_attr(vararr_xml, "name");
			fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", streamname, optname);
		}
	}

	for (var_xml = ezxml_child(varstruct_xml, "var"); var_xml; var_xml = var_xml->next) {
		optstream = ezxml_attr(var_xml, "streams");
		if (optstream != NULL && strstr(optstream, streamname) != NULL) {
			optname = ezxml_attr(var_xml, "name");
			fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", streamname, optname);
		}
	}

	for (substruct_xml = ezxml_child(varstruct_xml, "var_struct"); substruct_xml; substruct_xml = substruct_xml->next){
		generate_immutable_struct_contents(fd, streamname, substruct_xml);
	}

	return 0;
}/*}}}*/


/*********************************************************************************
 *
 *  Function: generate_immutable_streams
 *
 *  Generates the Fortran include file 'setup_immutable_streams.inc' that contains
 *  the subroutine mpas_generate_immutable_streams() responsible for making calls
 *  to the stream manager to define all immutable streams.
 *  The mpas_generate_immutable_streams() routine should be called after blocks
 *  have been allocated in the framework and after the stream manager has been
 *  initialized, but before any calls to generate mutable streams are made.
 *
 *********************************************************************************/
int generate_immutable_streams(ezxml_t registry){/*{{{*/
	ezxml_t streams_xml, stream_xml, var_xml, vararr_xml, varstruct_xml;
	ezxml_t substream_xml, matchstreams_xml, matchstream_xml;

	const char *optname, *opttype, *optvarname, *optstream, *optfilename, *optinterval_in, *optinterval_out, *optimmutable;
	const char *optstructname, *optsubstreamname, *optmatchstreamname, *optmatchimmutable;
	const char *corename;
	FILE *fd;

	char core_string[1024];

	corename = ezxml_attr(registry, "core");

	sprintf(core_string, "_%s_", corename);

	// For now, don't include core name in subroutines
	sprintf(core_string, "_");

	fd = fopen("setup_immutable_streams.inc", "w+");

	fprintf(stderr, "---- GENERATING IMMUTABLE STREAMS ----\n");

	fortprintf(fd, "subroutine mpas%ssetup_immutable_streams(manager)\n\n", core_string);
	fortprintf(fd, "   use MPAS_stream_manager, only : MPAS_streamManager_type, MPAS_STREAM_INPUT_OUTPUT, MPAS_STREAM_INPUT, &\n");
	fortprintf(fd, "                                   MPAS_STREAM_OUTPUT, MPAS_STREAM_NONE, MPAS_STREAM_PROPERTY_IMMUTABLE, &\n");
	fortprintf(fd, "                                   MPAS_stream_mgr_create_stream, MPAS_stream_mgr_add_field, MPAS_stream_mgr_set_property\n\n");
	fortprintf(fd, "   implicit none\n\n");
	fortprintf(fd, "   type (MPAS_streamManager_type), pointer :: manager\n\n");
	fortprintf(fd, "   integer :: ierr\n\n");

	for (streams_xml = ezxml_child(registry, "streams"); streams_xml; streams_xml = streams_xml->next) {
		for (stream_xml = ezxml_child(streams_xml, "stream"); stream_xml; stream_xml = stream_xml->next) {

			optimmutable = ezxml_attr(stream_xml, "immutable");

			if (optimmutable != NULL && strcmp(optimmutable, "true") == 0) {

				optname = ezxml_attr(stream_xml, "name");
				opttype = ezxml_attr(stream_xml, "type");
				optfilename = ezxml_attr(stream_xml, "filename_template");

				/* create the stream */
				if (strstr(opttype, "input") != NULL && strstr(opttype, "output") != NULL)
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_INPUT_OUTPUT, \'%s\', ierr=ierr)\n", optname, optfilename);
				else if (strstr(opttype, "input") != NULL)
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_INPUT, \'%s\', ierr=ierr)\n", optname, optfilename);
				else if (strstr(opttype, "output") != NULL)
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_OUTPUT, \'%s\', ierr=ierr)\n", optname, optfilename);
				else
					fortprintf(fd, "   call MPAS_stream_mgr_create_stream(manager, \'%s\', MPAS_STREAM_NONE, \'%s\', ierr=ierr)\n", optname, optfilename);

				/* Loop over streams listed within the stream (only use immutable streams) */
				for (substream_xml = ezxml_child(stream_xml, "stream"); substream_xml; substream_xml = ezxml_next(substream_xml)) {
					optsubstreamname = ezxml_attr(substream_xml, "name");

					/* Find stream definition with matching name */
					for (matchstreams_xml = ezxml_child(registry, "streams"); matchstreams_xml; matchstreams_xml = matchstreams_xml->next){
						for (matchstream_xml = ezxml_child(matchstreams_xml, "stream"); matchstream_xml; matchstream_xml = matchstream_xml->next){
							optmatchstreamname = ezxml_attr(matchstream_xml, "name");
							optmatchimmutable = ezxml_attr(matchstream_xml, "immutable");

							if (optmatchstreamname != NULL && strcmp(optmatchstreamname, optsubstreamname) == 0){
								if (optmatchimmutable != NULL && strcmp(optmatchimmutable, "true") == 0) {
									/* Loop over fields listed within the stream */
									for (var_xml = ezxml_child(matchstream_xml, "var"); var_xml; var_xml = var_xml->next) {
										optvarname = ezxml_attr(var_xml, "name");
										fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
									}

									/* Loop over arrays of fields listed within the stream */
									for (vararr_xml = ezxml_child(matchstream_xml, "var_array"); vararr_xml; vararr_xml = vararr_xml->next) {
										optvarname = ezxml_attr(vararr_xml, "name");
										fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
									}

									/* Loop over var structs listed within the stream */
									for (varstruct_xml = ezxml_child(matchstream_xml, "var_struct"); varstruct_xml; varstruct_xml = varstruct_xml->next) {
										optstructname = ezxml_attr(varstruct_xml, "name");
										fortprintf(fd, "   call MPAS_stream_mgr_add_pool(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optstructname);
									}

								} else {
									printf("ERROR: Immutable streams cannot contain mutable streams within them.\n");	
									printf("ERROR:     Immutable stream \'%s\' contains a mutable stream \'%s\'.\n", optname, optsubstreamname);
									return 1;
								}
							}
						}
					}
				}

				/* Loop over var structs listed within the stream */
				for (varstruct_xml = ezxml_child(stream_xml, "var_struct"); varstruct_xml; varstruct_xml = ezxml_next(varstruct_xml)) {
					optstructname = ezxml_attr(varstruct_xml, "name");
					fortprintf(fd, "   call MPAS_stream_mgr_add_pool(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optstructname);
				}


				/* Loop over arrays of fields listed within the stream */
				for (vararr_xml = ezxml_child(stream_xml, "var_array"); vararr_xml; vararr_xml = ezxml_next(vararr_xml)) {
					optvarname = ezxml_attr(vararr_xml, "name");
					fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
				}

				/* Loop over fields listed within the stream */
				for (var_xml = ezxml_child(stream_xml, "var"); var_xml; var_xml = ezxml_next(var_xml)) {
					optvarname = ezxml_attr(var_xml, "name");
					fortprintf(fd, "   call MPAS_stream_mgr_add_field(manager, \'%s\', \'%s\', ierr=ierr)\n", optname, optvarname);
				}

				/* Loop over fields looking for any that belong to the stream */
				for (varstruct_xml = ezxml_child(registry, "var_struct"); varstruct_xml; varstruct_xml = ezxml_next(varstruct_xml)) {
					generate_immutable_struct_contents(fd, optname, varstruct_xml);
				}

				fortprintf(fd, "   call MPAS_stream_mgr_set_property(manager, \'%s\', MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)\n\n", optname);
			}

		}
	}

	fortprintf(fd, "end subroutine mpas%ssetup_immutable_streams\n", core_string);

	fclose(fd);

	return 0;
}/*}}}*/


int push_attributes(ezxml_t currentPosition){/*{{{*/
	ezxml_t child_xml, child_xml2;
	ezxml_t childStruct1, childStruct2, lastStruct;

	const char *name, *name2;
	const char *subname;
	const char *super_time_levs, *super_packages;
	const char *sub_time_levs, *sub_packages;

	int skip_struct;

	name = ezxml_attr(currentPosition, "name");

	// Iterate over var_arrays
	for(child_xml = ezxml_child(currentPosition, "var_array"); child_xml; child_xml = child_xml->next){
		super_time_levs = ezxml_attr(currentPosition, "time_levs");
		super_packages = ezxml_attr(currentPosition, "packages");
		subname = ezxml_attr(child_xml, "name");
		sub_time_levs = ezxml_attr(child_xml, "time_levs");
		sub_packages = ezxml_attr(child_xml, "packages");

		if(!sub_time_levs && super_time_levs){
			child_xml = ezxml_set_attr(child_xml, "time_levs", super_time_levs);
		}

		if(!sub_packages && super_packages){
			child_xml = ezxml_set_attr(child_xml, "packages", super_packages);
		}

		// Iterate over vars in var_array
		for(child_xml2 = ezxml_child(child_xml, "var"); child_xml2; child_xml2 = child_xml2->next){
			super_packages = ezxml_attr(child_xml, "packages");
			sub_packages = ezxml_attr(child_xml2, "packages");

			if(!sub_packages && super_packages){
				child_xml2 = ezxml_set_attr(child_xml2, "packages", super_packages);
			}
		}
	}

	// Iterate over vars
	for(child_xml = ezxml_child(currentPosition, "var"); child_xml; child_xml = child_xml->next){
		super_packages = ezxml_attr(currentPosition, "packages");
		super_time_levs = ezxml_attr(currentPosition, "time_levs");
		subname = ezxml_attr(child_xml, "name");
		sub_time_levs = ezxml_attr(child_xml, "time_levs");
		sub_packages = ezxml_attr(child_xml, "packages");

		if(!sub_time_levs && super_time_levs){
			child_xml = ezxml_set_attr(child_xml, "time_levs", super_time_levs);
		}

		if(!sub_packages && super_packages){
			child_xml = ezxml_set_attr(child_xml, "packages", super_packages);
		}
	}

	// Iterate over var structs
	for(child_xml = ezxml_child(currentPosition, "var_struct"); child_xml; child_xml = child_xml->next){
		super_packages = ezxml_attr(currentPosition, "packages");
		super_time_levs = ezxml_attr(currentPosition, "time_levs");

		subname = ezxml_attr(child_xml, "name");
		sub_time_levs = ezxml_attr(child_xml, "time_levs");
		sub_packages = ezxml_attr(child_xml, "packages");

		if(!sub_time_levs && super_time_levs){
			child_xml = ezxml_set_attr(child_xml, "time_levs", super_time_levs);
		}

		if(!sub_packages && super_packages){
			child_xml = ezxml_set_attr(child_xml, "packages", super_packages);
		}

		push_attributes(child_xml);
	}

	return 0;
}/*}}}*/


int merge_structs_and_var_arrays(ezxml_t currentPosition){/*{{{*/
	ezxml_t old_child, new_child;
	ezxml_t childStruct1, childStruct2, lastStruct;

	const char *name, *name2;
	const char *subname;

	int skip_struct;

	// Merge var_structs
	for(childStruct1 = ezxml_child(currentPosition, "var_struct"); childStruct1; childStruct1 = childStruct1->next){
		name = ezxml_attr(childStruct1, "name");

		skip_struct = 0;
		for(childStruct2 = ezxml_child(currentPosition, "var_struct"); childStruct2 != childStruct1 && childStruct2; childStruct2 = childStruct2->next){
			name2 = ezxml_attr(childStruct2, "name");

			if(strcmp(name, name2) == 0){
				skip_struct = 1;
			}
		}

		if(!skip_struct && childStruct1->next){
			lastStruct = childStruct1;
			for(childStruct2 = childStruct1->next; childStruct2; childStruct2 = childStruct2->next){
				name2 = ezxml_attr(childStruct2, "name");

				if(strcmp(name, name2) == 0){
					// Merge children into childStruct1, and "remove" childStruct2
					for(old_child = ezxml_child(childStruct2, "var"); old_child; old_child = old_child->next){
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					for(old_child = ezxml_child(childStruct2, "var_array"); old_child; old_child = old_child->next){
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					for(old_child = ezxml_child(childStruct2, "var_struct"); old_child; old_child = old_child->next){
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					// Remove childStruct2
					lastStruct->next = childStruct2->next;
					free(childStruct2);
					childStruct2 = lastStruct;
				} else {
					lastStruct = childStruct2;
				}
			}
		}
	}

	// Merge var_arrays
	for(childStruct1 = ezxml_child(currentPosition, "var_array"); childStruct1; childStruct1 = childStruct1->next){
		name = ezxml_attr(childStruct1, "name");

		skip_struct = 0;
		for(childStruct2 = ezxml_child(currentPosition, "var_array"); childStruct2 && childStruct2 != childStruct1; childStruct2 = childStruct2->next){
			name2 = ezxml_attr(childStruct2, "name");

			if(strcmp(name, name2) == 0){
				skip_struct = 1;
			}
		}

		if(!skip_struct && childStruct1->next){
			lastStruct = childStruct1;
			for(childStruct2 = childStruct1->next; childStruct2; childStruct2 = childStruct2->next){
				name2 = ezxml_attr(childStruct2, "name");

				if(strcmp(name, name2) == 0){
					// Merge var_arrays and remove childStruct2
					for(old_child = ezxml_child(childStruct2, "var"); old_child; old_child = old_child->next){
						new_child = ezxml_insert(old_child, childStruct1, strlen(childStruct1->txt));
					}

					lastStruct->next = childStruct2->next;
					free(childStruct2);
					childStruct2 = lastStruct;
				} else {
					lastStruct = childStruct2;
				}
			}
		}
	}

	for(childStruct1 = ezxml_child(currentPosition, "var_struct"); childStruct1; childStruct1 = childStruct1->next){
		merge_structs_and_var_arrays(childStruct1);
	}

	return 0;
}/*}}}*/


int merge_streams(ezxml_t registry){/*{{{*/
	ezxml_t old_child, new_child, tmp_child;
	ezxml_t childStream1, childStream2, lastStream;
	ezxml_t includeStream;

	ezxml_t streamsBlock, streamsBlock2;

	const char *name, *name2;
	const char *subname;

	int skip_stream;

	// First, merge all streams blocks. Regardless of nested stream names.
	streamsBlock = ezxml_child(registry, "streams");
	while(streamsBlock->next){
		for(childStream1 = ezxml_child(streamsBlock->next, "stream"); childStream1; childStream1){
			if (childStream1->next){
				lastStream = childStream1->next;
			} else {
				lastStream = NULL;
			}
			name = ezxml_attr(childStream1, "name");
			new_child = ezxml_insert(childStream1, streamsBlock, strlen(streamsBlock->txt));

			childStream1 = lastStream;

		}

		lastStream = streamsBlock->next;
		streamsBlock->next = streamsBlock->next->next;
		free(lastStream);
	}


	// Now, merge all streams with the same name, within streamsBlock
	streamsBlock = ezxml_child(registry, "streams");
	for(childStream1 = ezxml_child(streamsBlock, "stream"); childStream1; childStream1 = childStream1->next){
		name = ezxml_attr(childStream1, "name");

		skip_stream = 0;

		for(childStream2 = ezxml_child(streamsBlock, "stream"); childStream2 && childStream2 != childStream1; childStream2 = childStream2->next){
			name2 = ezxml_attr(childStream2, "name");

			if(strcmp(name, name2) == 0){
				skip_stream = 1;
			}
		}

		if(!skip_stream && childStream1->next){
			lastStream = childStream1;
			for(childStream2 = childStream1->next; childStream2; childStream2 = childStream2->next){
				name2 = ezxml_attr(childStream2, "name");

				if(strcmp(name, name2) == 0){
					// Merge child vars
					for(old_child = ezxml_child(childStream2, "var"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}

					for(old_child = ezxml_child(childStream2, "var_array"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}

					for(old_child = ezxml_child(childStream2, "var_struct"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}

					for(old_child = ezxml_child(childStream2, "stream"); old_child; old_child){
						if(old_child->next){
							tmp_child = old_child->next;
						} else {
							tmp_child = NULL;
						}
						new_child = ezxml_insert(old_child, childStream1, strlen(childStream1->txt));

						old_child = tmp_child;
					}




					lastStream->next = childStream2->next;
					free(childStream2);
					childStream2 = lastStream;
				} else {
					lastStream = childStream2;
				}
			}
		}
	}

	return 0;
}/*}}}*/


int parse_structs_from_registry(ezxml_t registry)/*{{{*/
{
	ezxml_t structs_xml, var_arr_xml, var_xml;
	ezxml_t packages_xml, package_xml;

	const char *corename, *packagename, *structname, *structpackages;
	FILE *fd;
	int err;

	char core_string[1024];
	char spacing[1024];
	char package_list[2048];

	int no_packages;

	char *string, *tofree, *token;

	corename = ezxml_attr(registry, "core");

	sprintf(core_string, "_%s_", corename);

	// For now, don't include core name in subroutines.
	sprintf(core_string, "_");

	fd = fopen("structs_and_variables.inc", "w+");

	for (structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next){
		err = parse_struct(fd, registry, structs_xml, 0, '\0', corename);
	}

	fortprintf(fd, "   subroutine mpas_generate_structs(block, structPool, dimensionPool, packagePool)\n");
	fortprintf(fd, "      type (block_type), pointer, intent(inout) :: block\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: structPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(inout) :: dimensionPool\n");
	fortprintf(fd, "      type (mpas_pool_type), intent(in) :: packagePool\n");

	fortprintf(fd, "\n");

	for (structs_xml = ezxml_child(registry, "var_struct"); structs_xml; structs_xml = structs_xml->next){
		structname = ezxml_attr(structs_xml, "name");

		fortprintf(fd, "      call mpas_generate%spool_%s(block, structPool, dimensionPool, packagePool)\n", core_string, structname);

		fortprintf(fd, "\n");
	}


	fortprintf(fd, "   end subroutine mpas_generate_structs\n");

	fclose(fd);

	return 0;
}/*}}}*/


