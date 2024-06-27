/*

    Copyright (C) 2019 Alois Schloegl <alois.schloegl@ist.ac.at>
    This file is part of the "BioSig for C/C++" repository
    (biosig4c++) at http://biosig.sf.net/

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.


References:
    https://stackoverflow.com/questions/44378764/hello-tensorflow-using-the-c-api
    https://stackoverflow.com/questions/41688217/how-to-load-a-graph-with-tensorflow-so-and-c-api-h-in-c-language
    https://tebesu.github.io/posts/Training-a-TensorFlow-graph-in-C++-API
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tensorflow/c/c_api.h>
#include "mex.h"
//#include "matrix.h"


#ifdef tmwtypes_h
  #if (MX_API_VER<=0x07020000)
    typedef int mwSize;
  #endif
#endif


TF_Buffer* read_file(const char* file);

void free_buffer(void* data, size_t length) {
        free(data);
}

TF_Buffer* read_file(const char* file) {
  FILE *f = fopen(file, "rb");
  if (f==NULL) return NULL;
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  fseek(f, 0, SEEK_SET);  //same as rewind(f);

  void* data = malloc(fsize);
  fread(data, fsize, 1, f);
  fclose(f);

  TF_Buffer* buf = TF_NewBuffer();
  buf->data      = data;
  buf->length    = fsize;
  buf->data_deallocator = free_buffer;
  return buf;
}


void deallocTFData(void* data, size_t len, void* arg) {
	return;
}; 

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[]         /* array of pointers to input arguments */
)

{
	const mxArray	*arg;
	TF_Buffer* graph_def = NULL;
	TF_Tensor * tensor   = NULL;
		
	if (nrhs<1) {
                mexPrintf("mexTF (mexTensorflow) is in a very experimental state.\n");
		mexPrintf("   Usage of mexTF:\n");
		mexPrintf("\tv = mexTF()\n\t\treturns tensorflow version\n");
		mexPrintf("\t[v, graph_def] = mexTF('graph_def')\n\t\treads graph definition file\n");
		mexPrintf("\t[v, graph_def2] = mexTF(graph_def)\n\t\treads graph definition\n");
		mexPrintf("\t[v, graph_def2, class] = mexTF(graph_def, data)\n\t\treads graph definition\n");
		mexPrintf("   Input:\n");
		mexPrintf("   Output:\nTensorflow version\n");
        }

        mexPrintf("%s line %d: %d %d\n",__FILE__,__LINE__,nrhs,nlhs);
	for (int k = 0; k < nrhs; k++) {	
		arg = prhs[k];
		mxClassID argtype = mxGetClassID(arg);
		
		if (mxIsEmpty(arg) && (k>0)) {
			mexPrintf("%s line %d\n",__FILE__,__LINE__);
		}

		else if ( mxIsChar(arg) && (k==0) )  {
			mexPrintf("%s line %d\n",__FILE__,__LINE__);
			char *tmp = mxArrayToString(arg);
			graph_def = read_file(tmp);
			mxFree(tmp);
		}

		else if ( ((argtype==mxINT8_CLASS) || (argtype==mxUINT8_CLASS)) && (k==0) )  {
			mexPrintf("%s line %d\n",__FILE__,__LINE__);
			if (!graph_def) {
				graph_def = TF_NewBuffer();
				graph_def->data   = mxGetData(arg);
				graph_def->length = mxGetNumberOfElements(arg);
				graph_def->data_deallocator = NULL;
			};
		}

		else if ( (k==1) && mxIsNumeric(arg) )  {
			mexPrintf("%s line %d\n",__FILE__,__LINE__);
			TF_DataType tf_type; 

			mxClassID typ = mxGetClassID(arg);
			switch (argtype) {
			case mxDOUBLE_CLASS: 
				tf_type = TF_DOUBLE;
				break;
			case mxSINGLE_CLASS: 
				tf_type = TF_FLOAT;
				break;

			case mxINT64_CLASS: 
				tf_type = TF_INT64;
				break;
			case mxINT32_CLASS: 
				tf_type = TF_INT32;
				break;
			case mxINT16_CLASS: 
				tf_type = TF_INT16;
				break;
			case mxINT8_CLASS: 
				tf_type = TF_INT8;
				break;

			case mxUINT64_CLASS: 
				tf_type = TF_UINT64;
				break;
			case mxUINT32_CLASS: 
				tf_type = TF_UINT32;
				break;
			case mxUINT16_CLASS: 
				tf_type = TF_UINT16;
				break;
			case mxUINT8_CLASS: 
				tf_type = TF_UINT8;
				break;

			default:
				mexPrintf("Error: data type %s of arg1 not supported\n",mxGetClassName(arg));
				return;
				;
			}

			int ndims = mxGetNumberOfDimensions(arg);
			int64_t *dims = calloc(ndims, sizeof(int64_t));
			for (int k=0; k < ndims; k++) {
				dims[k] = *(mxGetDimensions(arg) + k);
				mexPrintf("%s line %d: dim[%d]= %d \n", __FILE__, __LINE__, k, dims[k]); 				
			}

			mexPrintf("%s line %d: going to converted to tensor [%d,%d,%d] \n", __FILE__, __LINE__, ndims, mxGetNumberOfElements(arg), TF_DataTypeSize(tf_type));

			tensor = TF_NewTensor( tf_type, dims, ndims, (void*)mxGetData(arg), mxGetNumberOfElements(arg) * TF_DataTypeSize(tf_type), &deallocTFData, NULL);

			mexPrintf("%s line %d: input converted to tensor %p\n", __FILE__, __LINE__, tensor);
			mexPrintf("%s line %d: input converted to tensor %d %d %d %d \n", __FILE__, __LINE__, TF_NumDims(tensor), TF_TensorByteSize(tensor), TF_Dim(tensor, 0), TF_Dim(tensor, 1));

			free(dims);
		}
        }

	mexPrintf("%s line %d\n",__FILE__,__LINE__);
        plhs[0] = mxCreateString(TF_Version());
	if ( (nlhs > 1) && graph_def ) {
		mexPrintf("%s line %d\n",__FILE__,__LINE__);
		const int ndim = 2;
		mwSize dims[ndim];
		dims[0] = 1;
		dims[1] = graph_def->length;
	        plhs[1] = mxCreateNumericArray(ndim, dims,  mxUINT8_CLASS, mxREAL);
		void *p = mxMalloc(dims[1]);
		memcpy(p, graph_def->data, dims[1]);
		mxSetData(plhs[1], p);
	}
	mexPrintf("%s line %d\n",__FILE__,__LINE__);

	/***********************************************
		load graph
	 ***********************************************/
	// Graph definition from unzipped https://storage.googleapis.com/download.tensorflow.org/models/inception5h.zip
	// which is used in the Go, Java and Android examples
	//   TF_Buffer* graph_def = read_file("inception5h/tensorflow_inception_graph.pb");
	TF_Graph* graph = TF_NewGraph();

	// Import graph_def into graph
	TF_Status* status = TF_NewStatus();
	TF_ImportGraphDefOptions* opts = TF_NewImportGraphDefOptions();
	TF_GraphImportGraphDef(graph, graph_def, opts, status);
	TF_DeleteImportGraphDefOptions(opts);
	TF_DeleteBuffer(graph_def);

	if (TF_GetCode(status) != TF_OK) {
		fprintf(stderr, "ERROR: Unable to import graph <%s>\n", TF_Message(status));
		TF_DeleteStatus(status);
		return;
	}
	fprintf(stdout, "Successfully imported graph\n");


if (tensor==NULL) {
	// Use the graph
	TF_DeleteGraph(graph);
	return;
}
	
	/***********************************************
		run session
	 ***********************************************/
	TF_SessionOptions * options = TF_NewSessionOptions();
	TF_Session * session = TF_NewSession( graph, options, status );

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	char hello[] = "Hello TensorFlow!";
//	if (tensor==NULL) tensor = TF_AllocateTensor( TF_STRING, 0, 0, 8 + TF_StringEncodedSize( strlen( hello ) ) );

	TF_Tensor * tensorOutput;

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	TF_OperationDescription * operationDescription = TF_NewOperation( graph, "Const", "hello" );

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	TF_Operation * operation;
	struct TF_Output output;

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

//	TF_StringEncode( hello, strlen( hello ), 8 + ( char * ) TF_TensorData( tensor ), TF_StringEncodedSize( strlen( hello ) ), status );
//	memset( TF_TensorData( tensor ), 0, 8 );

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	TF_SetAttrTensor( operationDescription, "value", tensor, status );

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	TF_SetAttrType( operationDescription, "dtype", TF_TensorType( tensor ) );

mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	operation = TF_FinishOperation( operationDescription, status );


mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	output.oper = operation;
	output.index = 0;

	TF_SessionRun( session, 0,
                 0, 0, 0,  // Inputs
                 &output, &tensorOutput, 1,  // Outputs
                 &operation, 1,  // Operations
                 0, status );


mexPrintf("%s line %d: %s\n",__FILE__,__LINE__, TF_Message(status));

	printf( "status code: %i\n", TF_GetCode( status ) );
	printf( "%s\n", ( ( char * ) TF_TensorData( tensorOutput ) ) + 9 );

	TF_CloseSession( session, status );
	TF_DeleteSession( session, status );
	TF_DeleteStatus( status );
	TF_DeleteSessionOptions( options );
	TF_DeleteGraph(graph);

}


