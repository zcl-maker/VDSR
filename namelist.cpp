#include "vdsr.h"
#include "Logging.h"

char NList::start_tag = '&';
char NList::close_tag = '/';
char NList::rem_tag = '#';

// --------------- Int   Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, int* ptr, int v)
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- Long   Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, long* ptr, long v)
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- Float Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, float* ptr, float v) 
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- Double Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, double* ptr, double v) 
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- String Entry ----------------------------------------

NList_Entry* NList::AddEntry(char* n, char* ptr, char* v) 
{
   NList* l = this;
   NList_Entry* e = new NList_Entry(l, n, ptr, v);
   e->next=head; head = e;
   return e;
}

// --------------- NList::NList ----------------------------------------

NList::NList() { // This constructor has been made for Controls::C_f_SaveTime.
//   NList::start_tag = '&';
//   NList::close_tag = '/';
//   NList::rem_tag = '#';

   tail=head=NULL;
   bsize=BUF_SIZE;
   err_msg[0] = "FILE_NOT_OPEN";
   err_msg[1] = "NLIST_NOT_FOUND";
   err_msg[2] = "WRONG_MEMBER";
   err_msg[3] = "CONVERSION_ERROR";
   err_msg[4] = "PREMATURE_END_OF_FILE";
}


// --------------- NList::NList ----------------------------------------

NList::NList(const char* name) {
   list_name = new char[strlen(name)+2];
   sprintf(list_name,"%s",name);
   tail=head=NULL; 
   bsize=BUF_SIZE;
   err_msg[0] = "FILE_NOT_OPEN";
   err_msg[1] = "NLIST_NOT_FOUND";
   err_msg[2] = "WRONG_MEMBER";
   err_msg[3] = "CONVERSION_ERROR";
   err_msg[4] = "PREMATURE_END_OF_FILE";
}

// --------------- Read NList ----------------------------------------

int NList::read(FILE* f) 
{
   char *stmp = NULL; 
   int found = -1;
  //   cout << "Looking for Namelist:" << list_name <<endl; getchar();

   if (f == NULL) error(1);
   while (fgets(buf,bsize,f)) {


      if (buf[0] == rem_tag) continue;
      if (stmp=strchr(buf,start_tag)) {
         stmp = strtok(stmp+1," \t\n");
         if (strcmp(stmp,list_name)) continue;
         break; 
      };
   };

   if (stmp==NULL) error(2);
   while (stmp=fgets(buf,bsize,f)) {
      if (buf[0] == rem_tag) continue;
      if (strchr(buf,close_tag)) {
         //      cout << "Namelist:" << list_name << " read" <<endl;
         return 0;
      }
      if ((stmp = strtok(stmp," =\t\n")) == NULL) continue;
      if (strcmp(stmp,"\n")==0) continue;
      NList_Entry* e = head;
      while (e) {
         if (strcmp(stmp,e->name)) { e=e->next; continue; }
         if ((stmp = strtok(NULL," =\t\n")) == NULL) continue;
         switch (e->t) {
      case 'i':  if (sscanf(stmp,"%d",(int*)e->ptr)==0) error(4); break;
      case 'l':  if (sscanf(stmp,"%d",(long*)e->ptr)==0) error(4); break;
      case 'f':  if (sscanf(stmp,"%g",(float*)e->ptr)==0) error(4); break;
      case 'd':  
         { float dtmp; double* dptr;
         if (sscanf(stmp,"%g",&dtmp)==0) error(4);
         dptr=(double*)e->ptr; *dptr = dtmp; break;
         }
      case 's':  if (sscanf(stmp,e->fmt,(char*)e->ptr)==0) error(4); 
         break;
      default: 
         ;}
         break;
      }
      if (stmp==NULL) continue;
      if (e==NULL) continue; //<SergK>{ cout << stmp << "\n"; error(3); }
   }

   error(5);
   return -5;
}

// --------------- Write NList ----------------------------------------

int NList::write(FILE* f) 
{
   if (f == NULL) error(1);
   fprintf(f,"%c%s \n",start_tag,list_name);
   NList_Entry* e = head;

   while (e) {
      fprintf(f,"%s = ",e->name);
      switch (e->t) {
      case 'i':  fprintf(f,"%d \n",*((int*)e->ptr)); break;
      case 'l':  fprintf(f,"%d \n",*((long*)e->ptr)); break;
      case 'f':  fprintf(f,"%g \n",*((float*)e->ptr)); break;
      case 'd':  fprintf(f,"%g \n",*((double*)e->ptr)); break;
      case 's':  
         fprintf(f,e->fmt,(char*)e->ptr);
         fprintf(f," \n");
         break;
      default: 
         ;
      }
      e=e->next;
   }

   fprintf(f,"%c \n",close_tag);

   return 0;
}


// --------------- Read NList ----------------------------------------

int NList::sread(char* f) 
{
   char *stmp = NULL; 
   int found = -1;
   //  cout << "Looking for Namelist:" << list_name <<endl;

   stmp = f;
   while( stmp=strchr(stmp,start_tag) )
   {
      sscanf(stmp+1,"%s",buf);

      if( !(strcmp(buf,list_name)) )
         break;

      stmp = stmp + 1;
   };

   f = stmp;
   if( stmp==NULL) error(2);

   stmp=strchr(f,'\n');
   f = stmp+1;

   while ( stmp=strchr(f,'\n') ) {
      int len = stmp-f;
      memcpy(buf,f,len);
      //<SergK/>
      buf[len] = '\0';
      //</SergK>
      f = stmp+1;

      if (buf[0] == rem_tag) continue;
      if (strchr(buf,close_tag))
      {
         //      cout << "Namelist:" << list_name << " read" <<endl;
         return 0;
      }
      if ((stmp = strtok(buf," =\t\n")) == NULL) continue;
      if (strcmp(stmp,"\n")==0) continue;
      NList_Entry* e = head;
      while (e)
      {
         //<SergK/>
         if( stmp==NULL)
         {
	    cout << "WARNING: NList::sread: stmp == NULL, e->name = " << e->name << endl;
            break;
         }
         //</SergK>         
         if (strcmp(stmp,e->name)) { e=e->next; continue; }
         if ((stmp = strtok(NULL," =\t\n")) == NULL) continue;
         switch (e->t) {
      case 'i':  if (sscanf(stmp,"%d",(int*)e->ptr)==0) error(4); break;
      case 'l':  if (sscanf(stmp,"%d",(long*)e->ptr)==0) error(4); break;
      case 'f':  if (sscanf(stmp,"%g",(float*)e->ptr)==0) error(4); break;
      case 'd':  
         { float dtmp; double* dptr;
         if (sscanf(stmp,"%g",&dtmp)==0) error(4);
         dptr=(double*)e->ptr; *dptr = dtmp; break;
         }
      case 's':  if (sscanf(stmp,e->fmt,(char*)e->ptr)==0) error(4); 
         break;
      default: 
         ;}
         break;
      }
      if (stmp==NULL) continue;
      if (e==NULL) { cout << stmp << "\n"; error(3); }
   }

   error(5);
   return -5;
}

// --------------- Write NList in a string ----------------------------------------

int NList::swrite(char* f) 
{
   int i=0;

   if (f == NULL) error(1);
   i += sprintf(f,"%c%s \n",start_tag,list_name);
   NList_Entry* e = head;

   while (e)
   {
      i += sprintf(f+i,"%s = ",e->name);
      switch (e->t) {
      case 'i':  i += sprintf(f+i,"%d \n",*((int*)e->ptr)); break;
      case 'l':  i += sprintf(f+i,"%d \n",*((long*)e->ptr)); break;
      case 'f':  i += sprintf(f+i,"%g \n",*((float*)e->ptr)); break;
      case 'd':  i += sprintf(f+i,"%g \n",*((double*)e->ptr)); break;
      case 's':  
         i += sprintf(f+i,e->fmt,(char*)e->ptr);
         i += sprintf(f+i," \n");
         break;
      default: 
         ;
      }
      e=e->next;
   }

   sprintf(f+i,"%c \n",close_tag);

   return 0;
}

// --------------- NList Error Handler ------------------------------------

void NList::error(int ierr) 
{
   cout << "NList "<< list_name<<" error: "<<err_msg[ierr-1] <<endl;
   exit(ierr);
}

// --------------- NList Packer ------------------------------------
void NList::pack_nls(CBuffer *buf)
{
   int tmp_size=0;
   NList_Entry* tmp = head;
   while (tmp)
   {
      tmp_size = tmp->size;
      buf->npush(tmp->ptr,tmp_size,1);
      tmp = tmp->next;
   }
}

// --------------- NList UnPacker ------------------------------------
void NList::unpack_nls(CBuffer *buf)
{
   NList_Entry* tmp = head;

   while (tmp) {
      buf->npop(tmp->ptr,tmp->size,1);
      tmp = tmp->next;
   }
   buf->reset();
}

// --------------- NList Constructors ------------------------------------

NList_Entry::NList_Entry(NList* l, char* n, char* p, char* v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   fmt = v;
   t = 's';
   int m_size = 0;
   if (sscanf(v+1,"%d",&size)) {
      m_size = size;
   } else {
      size = 1;
   };
}

NList_Entry::NList_Entry(NList* l, char* n, int* p, int v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = v;
   fmt = "%d";
   t = 'i';
   size = sizeof(int);
}

NList_Entry::NList_Entry(NList* l, char* n, long* p, long v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = v;
   fmt = "%d";
   t = 'l';
   size = sizeof(long);
}

NList_Entry::NList_Entry(NList* l, char* n, float* p, float v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = (float)v;
   fmt = "%g";
   t = 'f';
   size = sizeof(float);
}

NList_Entry::NList_Entry(NList* l, char* n, double* p, double v)
{
   name = new char[strlen(n)+2];
   sprintf(name,"%s",n);
   list = l;
   ptr = (void*)p;
   *p = (double)v;
   fmt = "%g";
   t = 'd';
   size = sizeof(double);
}
