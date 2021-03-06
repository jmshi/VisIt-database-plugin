/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            avtATHENAFileFormat.C                           //
// ************************************************************************* //

#include <avtATHENAFileFormat.h>

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>

#include <DebugStream.h>
#include <InvalidVariableException.h>
#include <InvalidFilesException.h>
#include <UnexpectedValueException.h>



using     std::string;


// ****************************************************************************
//  Method: avtATHENAFileFormat constructor
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

avtATHENAFileFormat::avtATHENAFileFormat(const char *filename)
    : avtSTMDFileFormat(&filename, 1)
{
// 1. open file
  OpenFile(filename);
  if(!fopened) {
	std::cout << "file handler is not returned properly in constructor" << std::endl;
  }
// 2. read the header
  fread(&time,sizeof(float),1,fh);
  std::cout << "time = " << time << std::endl;
  for (int n=0;n<8;n++) {
	fread(&ndata[n],sizeof(int),1,fh);
	//std::cout << "ndata[" << n <<"] = " << ndata[n] << std::endl;
  }
// 3. allocate memory and read for the coord
  x1 = new float[ndata[1]];
  fread(&x1[0],sizeof(float),ndata[1],fh);
  x2 = new float[ndata[2]];
  fread(&x2[0],sizeof(float),ndata[2],fh);
  x3 = new float[ndata[3]];
  fread(&x3[0],sizeof(float),ndata[3],fh);
  //for (int i=0;i<ndata[3];i++) {
  //  std::cout << "x3["<< i << "] = " << x3[i] << std::endl;
  //}

// 4. allocate memory and read for the data var

  int asize = ndata[6]*ndata[5]*ndata[4];
  var = new float[asize*ndata[7]];

  int p = 0;
  for (int n=0;n<ndata[7];n++) {
  for (int k=0;k<ndata[6];k++) {
  for (int j=0;j<ndata[5];j++) {
  for (int i=0;i<ndata[4];i++) {
	//fread(&(var(k,j,i)),sizeof(float),1,fh);
	// std::cout << std::setprecision(10) << var(k,j,i) << std::endl;
	fread(&(var[p++]),sizeof(float),1,fh);
	if(n==0 && k==50 && j==60) {
	  std::cout << std::setprecision(10) << var[p-1] << std::endl;
	}
  }}}}

  fclose(fh);

}


// ****************************************************************************
//  Method: avtATHENAFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

void
avtATHENAFileFormat::FreeUpResources(void)
{
}


// ****************************************************************************
//  Method: avtATHENAFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

void
avtATHENAFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
	std::cout << "file opened = " << fopened << std::endl;
    //
    // CODE TO ADD A MESH
	//
    string meshname = "Mesh";
    //
    // AVT_RECTILINEAR_MESH, AVT_CURVILINEAR_MESH, AVT_UNSTRUCTURED_MESH,
    // AVT_POINT_MESH, AVT_SURFACE_MESH, AVT_UNKNOWN_MESH
    avtMeshType mt = AVT_RECTILINEAR_MESH;
    //
    int nblocks = ndata[0];
    int block_origin = 0;
    int spatial_dimension = 3;
    int topological_dimension = 3;
    double *extents = NULL;
    //
    // Here's the call that tells the meta-data object that we have a mesh:
    //
    AddMeshToMetaData(md, meshname, mt, extents, nblocks, block_origin,
                       spatial_dimension, topological_dimension);


    //
    // CODE TO ADD A SCALAR VARIABLE
    //
    string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    string varname = "density";
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    avtCentering cent = AVT_ZONECENT;
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
    //
    varname = "vel1";
    AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
    varname = "vel2";
    AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
    varname = "vel3";
    AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
	if (ndata[7] > 5) {
      varname = "B1";
      AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
      varname = "B2";
      AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
      varname = "B3";
      AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
	}

    //
    // CODE TO ADD A VECTOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int vector_dim = 2;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddVectorVarToMetaData(md, varname, mesh_for_this_var, cent,vector_dim);
    //

    //
    // CODE TO ADD A TENSOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int tensor_dim = 9;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddTensorVarToMetaData(md, varname, mesh_for_this_var, cent,tensor_dim);
    //

    //
    // CODE TO ADD A MATERIAL
    //
    // string mesh_for_mat = meshname; // ??? -- could be multiple meshes
    // string matname = ...
    // int nmats = ...;
    // vector<string> mnames;
    // for (int i = 0 ; i < nmats ; i++)
    // {
    //     char str[32];
    //     sprintf(str, "mat%d", i);
    //     -- or --
    //     strcpy(str, "Aluminum");
    //     mnames.push_back(str);
    // }
    //
    // Here's the call that tells the meta-data object that we have a mat:
    //
    // AddMaterialToMetaData(md, matname, mesh_for_mat, nmats, mnames);
    //
    //
    // Here's the way to add expressions:
    //Expression momentum_expr;
    //momentum_expr.SetName("momentum");
    //momentum_expr.SetDefinition("{u, v}");
    //momentum_expr.SetType(Expression::VectorMeshVar);
    //md->AddExpression(&momentum_expr);
    //Expression KineticEnergy_expr;
    //KineticEnergy_expr.SetName("KineticEnergy");
    //KineticEnergy_expr.SetDefinition("0.5*(momentum*momentum)/(rho*rho)");
    //KineticEnergy_expr.SetType(Expression::ScalarMeshVar);
    //md->AddExpression(&KineticEnergy_expr);
    //
}


// ****************************************************************************
//  Method: avtATHENAFileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      domain      The index of the domain.  If there are NDomains, this
//                  value is guaranteed to be between 0 and NDomains-1,
//                  regardless of block origin.
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

vtkDataSet *
avtATHENAFileFormat::GetMesh(int domain, const char *meshname)
{
    //YOU MUST IMPLEMENT THIS
    int ndims = 3;
    int dims[3] = {1,1,1};
    vtkFloatArray *coords[3] = {0,0,0};

    // Read the ndims and number of X,Y,Z nodes from file.
    ndims = 3;
    dims[0] = ndata[1];
    dims[1] = ndata[2];
    dims[2] = ndata[3];

    // Read the X coordinates from the file.
    coords[0] = vtkFloatArray::New();
    coords[0]->SetNumberOfTuples(dims[0]);
    float *xarray = (float *)coords[0]->GetVoidPointer(0);
    //memcpy(xarray,x1,dims[0]);
	for (int i=0;i<dims[0];i++) {
	  *xarray = x1[i];
	  std::cout << "coords[0][" << i << "] = " << *(xarray) << std::endl;
	  xarray++;
	}

    // Read the Y coordinates from the file.
    coords[1] = vtkFloatArray::New();
    coords[1]->SetNumberOfTuples(dims[1]);
    float *yarray = (float *)coords[1]->GetVoidPointer(0);
    //memcpy(yarray,x2,dims[1]);
	for (int i=0;i<dims[1];i++) {
	  *yarray = x2[i];
	  yarray++;
	}

    // Read the Z coordinates from the file.
    coords[2] = vtkFloatArray::New();
    coords[2]->SetNumberOfTuples(dims[2]);
    float *zarray = (float *)coords[2]->GetVoidPointer(0);
    //memcpy(zarray,x3,dims[2]);
	for (int i=0;i<dims[2];i++) {
	  *zarray = x3[i];
	  zarray++;
	}

    //
    // Create the vtkRectilinearGrid object and set its dimensions
    // and coordinates.
    //
    vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
    rgrid->SetDimensions(dims);
    rgrid->SetXCoordinates(coords[0]);
    coords[0]->Delete();
    rgrid->SetYCoordinates(coords[1]);
    coords[1]->Delete();
    rgrid->SetZCoordinates(coords[2]);
    coords[2]->Delete();

    return rgrid;
	return 0;
}


// ****************************************************************************
//  Method: avtATHENAFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

vtkDataArray *
avtATHENAFileFormat::GetVar(int domain, const char *varname)
{
    int nvals;
    // Read the number of values contained in the array
    // specified by varname.
    nvals = ndata[4]*ndata[5]*ndata[6];

    // Allocate the return vtkFloatArray object. Note that
    // you can use vtkFloatArray, vtkDoubleArray,
    // vtkUnsignedCharArray, vtkIntArray, etc.
    vtkFloatArray *arr = vtkFloatArray::New();
    arr->SetNumberOfTuples(nvals);
    float *data = (float *)arr->GetVoidPointer(0);
    //READ nvals FLOAT NUMBERS INTO THE data ARRAY.
	int offset;
    if(strcmp(varname, "density")==0) {
 	  offset = 0;
	  std::cout << "load density offset = " << offset << std::endl;
	}
    if(strcmp(varname, "vel1")==0)    offset = nvals;
    if(strcmp(varname, "vel2")==0)    offset = nvals*2;
    if(strcmp(varname, "vel3")==0)    offset = nvals*3;
    if(strcmp(varname, "B1")==0)    offset = nvals*(ndata[7]-3);
    if(strcmp(varname, "B2")==0)    offset = nvals*(ndata[7]-2);
    if(strcmp(varname, "B3")==0)    offset = nvals*(ndata[7]-1);

    for (int i =0;i<100;i++)
	  std::cout << var[i] << " ";

	int p=0;
	for(int k=0;k<ndata[6];k++){
	for(int j=0;j<ndata[5];j++){
	for(int i=0;i<ndata[4];i++){
	   //std::cout << "inside loop [k,j,i] = " << k << " " << j << " " << i << std::endl;
	   //std::cout << ndata[6] << ndata[5] << ndata[4] << std::endl;
	   //std::cout << "var[" << offset+p << "] = " << var[offset+p] << std:: endl;
	  *data = var[offset+p];
	  if(k==30 && j==30)
		std::cout << "var[30,30,"<< i << "]=" << var[offset+p] << std::endl;
	  p++;
	  data++;
	}}}

    return arr;
	return 0;

    //
    // If you do have a scalar variable, here is some code that may be helpful.
    //
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // rv->SetNumberOfTuples(ntuples);
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      rv->SetTuple1(i, VAL);  // you must determine value for ith entry.
    // }
    //
    // return rv;
    //
}


// ****************************************************************************
//  Method: avtATHENAFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

vtkDataArray *
avtATHENAFileFormat::GetVectorVar(int domain, const char *varname)
{
    //YOU MUST IMPLEMENT THIS
	return 0;
    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a vector variable, here is some code that may be helpful.
    //
    // int ncomps = YYY;  // This is the rank of the vector - typically 2 or 3.
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // int ucomps = (ncomps == 2 ? 3 : ncomps);
    // rv->SetNumberOfComponents(ucomps);
    // rv->SetNumberOfTuples(ntuples);
    // float *one_entry = new float[ucomps];
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      int j;
    //      for (j = 0 ; j < ncomps ; j++)
    //           one_entry[j] = ...
    //      for (j = ncomps ; j < ucomps ; j++)
    //           one_entry[j] = 0.;
    //      rv->SetTuple(i, one_entry);
    // }
    //
    // delete [] one_entry;
    // return rv;
    //
}

// ****************************************************************************
//  Method: avtATHENAFileFormat::OpenFile
//
//  Purpose: open file and return a handler
//
//  Arguments:
//
//  Programmer: Ji-Ming Shi -- generated by xml2avt
//  Creation:   Fri Sep 9 16:01:29 PST 2016
//
// ****************************************************************************

void
avtATHENAFileFormat::OpenFile(const char *filename)
{
  if(filename==NULL) {
    // no input file is given
	debug4 << "### FATAL ERROR in OpenFile"
    << "No input file is specified." << std::endl;

	EXCEPTION1(InvalidFilesException, filename);
  }

  if((fh=fopen(filename,"rb"))==NULL) {
	debug4 << "### FATAL ERROR in OpenFile"
    << "unable to open input file." << std::endl;
	EXCEPTION1(InvalidFilesException, filename);
  }

  fopened = true;

}
