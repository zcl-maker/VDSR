
/* HDF5 Library                                                       */ 
#include "myhdfshell.h"
#include "vdsr.h"
int WriteHDF5RecordDouble(hid_t file, char* dataname, int nfields, double *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dimsfi[0] = nfields;
   dataspace = H5Screate_simple(rank1, dimsfi, NULL); 
   assert (dataspace >= 0);
   dataset = H5Dcreate(file, dataname, H5T_IEEE_F64LE, dataspace, H5P_DEFAULT);
   assert (dataset >= 0);
   status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, source);
   assert (status >= 0);
   status = H5Dclose(dataset);
   assert (status >= 0);
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}


int WriteHDF5RecordFloat(hid_t file, char* dataname, int nfields, float *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dimsfi[0] = nfields;
   dataspace = H5Screate_simple(rank1, dimsfi, NULL); 
   assert (dataspace >= 0);
   dataset = H5Dcreate(file, dataname, H5T_IEEE_F32LE, dataspace, H5P_DEFAULT);
   assert (dataset >= 0);
   status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, source);
   assert (status >= 0);
   status = H5Dclose(dataset);
   assert (status >= 0);
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}


int WriteHDF5RecordInt(hid_t file, char* dataname, int nfields, int *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dimsfi[0] = nfields;
   dataspace = H5Screate_simple(rank1, dimsfi, NULL); 
   assert (dataspace >= 0);
   dataset = H5Dcreate(file, dataname, H5T_STD_I32LE, dataspace, H5P_DEFAULT);
   assert (dataset >= 0);
   status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, source);
   assert (status >= 0);
   status = H5Dclose(dataset);
   assert (status >= 0);
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}

int WriteHDF5RecordLong(hid_t file, char* dataname, int nfields, long *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dimsfi[0] = nfields;
   dataspace = H5Screate_simple(rank1, dimsfi, NULL); 
   assert (dataspace >= 0);
   dataset = H5Dcreate(file, dataname, H5T_STD_I64LE, dataspace, H5P_DEFAULT);
   assert (dataset >= 0);
   status = H5Dwrite(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, source);
   assert (status >= 0);
   status = H5Dclose(dataset);
   assert (status >= 0);
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}

int WriteHDF5RecordString(hid_t file, char* dataname, int nfields, char *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dimsfi[0] = nfields;
   dataspace = H5Screate_simple(rank1, dimsfi, NULL); 
   assert (dataspace >= 0);
   dataset = H5Dcreate(file, dataname, H5T_STD_U8LE, dataspace, H5P_DEFAULT);
   assert (dataset >= 0);
   status = H5Dwrite(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, source);
   assert (status >= 0);
   status = H5Dclose(dataset);
   assert (status >= 0);
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}

int WriteHDF5Record(hid_t file, char* dataname, int nfields, double *source)
{
   return WriteHDF5RecordDouble(file, dataname, nfields, source);
}


int WriteHDF5Record(hid_t file, char* dataname, int nfields, float *source)
{
   return WriteHDF5RecordFloat(file, dataname, nfields, source);
}

int WriteHDF5Record(hid_t file, char* dataname, int nfields, int *source)
{
   return WriteHDF5RecordInt(file, dataname, nfields, source);
}

int WriteHDF5Record(hid_t file, char* dataname, int nfields, long *source)
{
   return WriteHDF5RecordLong(file, dataname, nfields, source);
}

int WriteHDF5Record(hid_t file, char* dataname, int nfields, char *source)
{
   return WriteHDF5RecordString(file, dataname, nfields, source);
}

int ReadHDF5Record(hid_t file, char* dataname, int nfields, float *source)
{
   return ReadHDF5RecordFloat(file, dataname, nfields, source);
}

int ReadHDF5Record(hid_t file, char* dataname, int nfields, int *source)
{
   return ReadHDF5RecordInt(file, dataname, nfields, source);
}

int ReadHDF5Record(hid_t file, char* dataname, int nfields, long *source)
{
   return ReadHDF5RecordLong(file, dataname, nfields, source);
}

int ReadHDF5Record(hid_t file, char* dataname, int nfields, char *source)
{
   return ReadHDF5RecordString(file, dataname, nfields, source);
}

int ReadHDF5RecordString(hid_t file, char* dataname, int nfields, char *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dataset = H5Dopen(file, dataname);
   assert (dataset >= 0);

   dataspace = H5Dget_space(dataset);   
   hid_t dummy_rank      = H5Sget_simple_extent_ndims(dataspace);
   assert(dummy_rank == rank1);
   status  = H5Sget_simple_extent_dims(dataspace, dimsfi, NULL);
   assert (status >= 0);
   assert(dimsfi[0] == nfields);
   status = H5Dread(dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, source);
   assert (status >= 0);
                                    // Close the dataspace                                            
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}

int ReadHDF5RecordInt(hid_t file, char* dataname, int nfields, int *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dataset = H5Dopen(file, dataname);
   assert (dataset >= 0);

   dataspace = H5Dget_space(dataset);   
   hid_t dummy_rank      = H5Sget_simple_extent_ndims(dataspace);
   assert(dummy_rank == rank1);
   status  = H5Sget_simple_extent_dims(dataspace, dimsfi, NULL);
   assert (status >= 0);
   assert(dimsfi[0] == nfields);
   status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, source);
   assert (status >= 0);
      // Close the dataspace                                            
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}

int ReadHDF5RecordLong(hid_t file, char* dataname, int nfields, long *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dataset = H5Dopen(file, dataname);
   assert (dataset >= 0);

   dataspace = H5Dget_space(dataset);   
   hid_t dummy_rank      = H5Sget_simple_extent_ndims(dataspace);
   assert(dummy_rank == rank1);
   status  = H5Sget_simple_extent_dims(dataspace, dimsfi, NULL);
   assert (status >= 0);
   assert(dimsfi[0] == nfields);
   status = H5Dread(dataset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, source);
   assert (status >= 0);
      // Close the dataspace                                            
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}

int ReadHDF5RecordFloat(hid_t file, char* dataname, int nfields, float *source)
{
   hid_t rank1 = 1;
   hsize_t dimsfi[3];
   hid_t       dataset;               /* File and dataset            */
   hid_t       dataspace;             /* Dataspace handles           */
   herr_t      status;                /* Error checking              */

   assert (file >= 0);

   dataset = H5Dopen(file, dataname);
   assert (dataset >= 0);

   dataspace = H5Dget_space(dataset);   
   hid_t dummy_rank      = H5Sget_simple_extent_ndims(dataspace);
   assert(dummy_rank == rank1);
   status  = H5Sget_simple_extent_dims(dataspace, dimsfi, NULL);
   assert (status >= 0);
   assert(dimsfi[0] == nfields);
   status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, source);
   assert (status >= 0);
      // Close the dataspace                                            
   status = H5Sclose(dataspace);
   assert (status >= 0);
   int istatus = status;
   return istatus;
}
