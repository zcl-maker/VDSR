/***************************************************************************
                          Min Chen copied from VLPL
                          Logging.h  -  description
                             -------------------
    begin                : Wed Jun 2 2004
    copyright            : (C) 2004 by Sergey Kiselev
    email                : sergk@thphy.uni-duesseldorf.de
 ***************************************************************************/
#ifndef H_LOGGING
#define H_LOGGING
#include "vdsr.h"
#include <iostream>
#include <string>

namespace Log
{
   class Logging
   {
      std::string _str;
      
      friend const char* ConstructorError(const std::string&);
      friend const char* ConstructorWarning(const std::string&);
      friend const char* Error(const std::string&);
      friend const char* Warning(const std::string&);
      friend const char* ConstructorErrorNoIni(const std::string&);
      friend const char* ConstructorErrorInputValues(const std::string&);
      friend const char* ConstructorErrorSetFunctionPreForParticle(const std::string&);
      friend const char* ConstructorErrorSetFunctionForParticle(const std::string&);
      friend const char* ConstructorErrorSetFunctionPostForParticle(const std::string&);
      friend const char* ConstructorErrorSetFunction(const std::string&);
      friend const char* ConstructorErrorSetFunctionPre(const std::string&);
      friend const char* ConstructorErrorInitializeAdoptedData(const std::string&);    

      public:
      //typedef const Log::Logging &(p_endll)(const Log::Logging &);
      typedef void(p_endl)();
      
      template <class Cl> inline const Logging& operator<<(Cl& cl) const
      {
         Domain::p_D->Getout_Flog() << cl;
         cout << cl;
         return *this;
      }
      template <class Cl> inline const Logging& operator<<(const Cl& cl) const
      {
         Domain::p_D->Getout_Flog() << cl;
         cout << cl;
         return *this;
      }      
      //template <p_endl> inline Logging& operator<< (p_endl& cl)
//      inline const Logging& operator<< (p_endll& cl) const
//      {
//         Domain::p_D->Getout_Flog() << std::endl;
//         std::cout << std::endl;
//         return *this;
//      }
      inline const Logging& operator<< (p_endl&) const
      {
         Domain::p_D->Getout_Flog() << std::endl;
         std::cout << std::endl;
         return *this;
      }
   };

//const Log::Logging& endl(const Log::Logging &);
void endl();
const char* ConstructorError(const std::string&);
const char* ConstructorWarning(const std::string&);
const char* Error(const std::string&);
const char* Warning(const std::string&);
const char* ConstructorErrorNoIni(const std::string&);
const char* ConstructorErrorInputValues(const std::string&);
const char* ConstructorErrorSetFunctionPreForParticle(const std::string&);
const char* ConstructorErrorSetFunctionForParticle(const std::string&);
const char* ConstructorErrorSetFunctionPostForParticle(const std::string&);
const char* ConstructorErrorSetFunction(const std::string&);
const char* ConstructorErrorSetFunctionPre(const std::string&);
const char* ConstructorErrorInitializeAdoptedData(const std::string&); 

}

extern Log::Logging logging;

#endif
