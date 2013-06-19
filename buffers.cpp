#include "buffers.h"

//---------------------- Class CBuffer::CBuffer -----------------------------
CBuffer::CBuffer(long size)
{
   buf = NULL;
   buf = new char[size];
   b_len = size;
   b_pos = 0;
   b_out = 0;
   b_CRC = 0;
}

//---CBuffer::------------->
bool CBuffer::Adjust(long size)
{
   while (b_len < size)
   {
      b_len += BUF_SIZE;
      
      char *btmp = new char[b_len];
      if(!btmp)
         return false;

      memcpy (btmp, buf, b_pos);
      delete[] buf;
      buf = btmp;
   }
   return true;
}

//---CBuffer::------------->
long CBuffer::FindCRC(long upto)
{
  if (upto<=0) upto = b_pos;
   b_CRC = 0;
   for (long i=0; i<upto; i++)
   {
     char c = buf[i];
     b_CRC += long(c);
   }
   return b_CRC;
}

//---CBuffer::------------->
#ifdef _LOAD_PARTICLES
void CBuffer::npush(void *ptr, long size, long nitems) throw(bad_alloc)
#else //_LOAD_PATICLES
void CBuffer::npush(void *ptr, long size, long nitems)
#endif//_LOAD_PATICLES
{
   long to_push = size*nitems;

#ifdef _LOAD_PARTICLES   
   if(!Adjust(b_pos+to_push))
      throw bad_alloc();
#else//_LOAD_PATICLES
   Adjust(b_pos+to_push);
#endif//_LOAD_PATICLES    

   if (ptr)
      memcpy (buf+b_pos, ptr, to_push);
   b_pos += to_push;
}

//---CBuffer::------------->
#ifdef _LOAD_PARTICLES
void CBuffer::npop(void *ptr, long size, long nitems) throw(CBufferError)
#else
void CBuffer::npop(void *ptr, long size, long nitems)
#endif//_LOAD_PATICLES
{
   long to_pop = size*nitems;

   if (ptr)
      memcpy (ptr, buf+b_out, to_pop);
   //  cout << "Extracting "<<int(*(buf+b_out))<<"\n";

   b_out += to_pop;
#ifdef _LOAD_PARTICLES   
   if (b_out == b_pos)
   {
      char message[128];
      sprintf(message," b_pos == b_out\n");
      std::string str_tmp = "CBuffer::npop";
      str_tmp += message;
      throw CBufferError(str_tmp);
   }
#else   
   if (b_out == b_pos)
   { 
      b_pos=b_out=0;
      return;
   }
#endif//_LOAD_PARTICLES   
   if (b_out > b_pos)
   { 
      cout << "Error CBuffer: b_out = "<< b_out<<" > b_pos = "<< b_pos<<"\n";
      exit(-3);
   }
   return;
}

//---CBuffer:: --------------------->
void CBuffer::copy(CBuffer* rbuf)
{
   reset();
   //  rbuf->dump();
   npush(rbuf->buf,rbuf->b_pos,1);
   rbuf->reset();
}

//---CBuffer::------------->
void CBuffer::print_state(void)
{
   cout << "CBuffer: len="<<b_len<<" pos="<<b_pos<<" out="<<b_out<<"\n";
}

//---CBuffer::------------->
void CBuffer::dump(void)
{
   long tmp = 0;
   cout <<"DUMP of CBuffer, pos = "<<b_pos<<"\n";
   for (tmp=0; tmp<b_pos; tmp+=4)
   {
      cout << "CBuffer["<<tmp<<"] ="<<int(*(buf+tmp))<<"\n";
   }
}

#undef BUF_SIZE

//-------------------- Class FCBuffer::FCBuffer -----------------------------
FCBuffer::FCBuffer(char *fname)
{
   filename = fname;
   file.open(fname, ios::in|ios::out);
   if (!file)
      cout << "ERROR FCBuffer:: Can not open file "<<fname<<endl;
   b_len = 0;
   b_pos = 0;
   b_out = 0;
}

//---FCBuffer::------------->
#ifdef _LOAD_PARTICLES
void FCBuffer::npush(void *ptr, long size, long nitems) throw(bad_alloc)
#else
void FCBuffer::npush(void *ptr, long size, long nitems)
#endif//_LOAD_PATICLES
{
   long to_push = size*nitems;
   if (file.good())
   {
      file.write((char*)ptr, to_push);
      b_len += to_push;
      b_pos += to_push;
   }

   else
      cout << "ERROR FCBuffer:: Can not WRITE "<<
      b_len+to_push<<" bytes in file "<<filename<<endl;
}

//---FCBuffer::------------->
#ifdef _LOAD_PARTICLES
void FCBuffer::npop(void *ptr, long size, long nitems) throw(CBufferError)
#else
void FCBuffer::npop(void *ptr, long size, long nitems)
#endif//_LOAD_PATICLES
{
   long to_pop = size*nitems;
   if (file.good())
   {
      file.read((char*)ptr, to_pop);
      b_out += to_pop;
   }
   else
      cout << "ERROR FCBuffer:: Can not READ "<<
      b_out+to_pop<<" bytes from file "<<filename<<endl;
   //  cout << "Extracting "<<int(*(buf+b_out))<<"\n";
   return;
}

//---FCBuffer:: --------------------->
void FCBuffer::reset(void)
{
   b_len = b_pos = b_out = 0;
   file.seekg(0);
   file.seekp(0);
}

//---FCBuffer::------------->
void FCBuffer::print_state(void)
{
   cout << "FCBuffer: filename="<<filename<<endl;
   cout << "FCBuffer: len="<<b_len<<" pos="<<b_pos<<" out="<<b_out<<endl;
}
