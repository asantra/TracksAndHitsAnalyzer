// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIarkasantradIarkadISasha_WorkdIlxmacrodIMHists_C_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/arkasantra/arka/Sasha_Work/lxmacro/./MHists.C"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *MHists_Dictionary();
   static void MHists_TClassManip(TClass*);
   static void *new_MHists(void *p = 0);
   static void *newArray_MHists(Long_t size, void *p);
   static void delete_MHists(void *p);
   static void deleteArray_MHists(void *p);
   static void destruct_MHists(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MHists*)
   {
      ::MHists *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MHists));
      static ::ROOT::TGenericClassInfo 
         instance("MHists", "MHists.h", 23,
                  typeid(::MHists), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MHists_Dictionary, isa_proxy, 4,
                  sizeof(::MHists) );
      instance.SetNew(&new_MHists);
      instance.SetNewArray(&newArray_MHists);
      instance.SetDelete(&delete_MHists);
      instance.SetDeleteArray(&deleteArray_MHists);
      instance.SetDestructor(&destruct_MHists);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MHists*)
   {
      return GenerateInitInstanceLocal((::MHists*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MHists*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MHists_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::MHists*)0x0)->GetClass();
      MHists_TClassManip(theClass);
   return theClass;
   }

   static void MHists_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/arkasantra/arka/Sasha_Work/lxmacro/./MHists.h");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_MHists(void *p) {
      return  p ? new(p) ::MHists : new ::MHists;
   }
   static void *newArray_MHists(Long_t nElements, void *p) {
      return p ? new(p) ::MHists[nElements] : new ::MHists[nElements];
   }
   // Wrapper around operator delete
   static void delete_MHists(void *p) {
      delete ((::MHists*)p);
   }
   static void deleteArray_MHists(void *p) {
      delete [] ((::MHists*)p);
   }
   static void destruct_MHists(void *p) {
      typedef ::MHists current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MHists

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 469,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::__1::vector<int, std::__1::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlETLegendmUgR_Dictionary();
   static void vectorlETLegendmUgR_TClassManip(TClass*);
   static void *new_vectorlETLegendmUgR(void *p = 0);
   static void *newArray_vectorlETLegendmUgR(Long_t size, void *p);
   static void delete_vectorlETLegendmUgR(void *p);
   static void deleteArray_vectorlETLegendmUgR(void *p);
   static void destruct_vectorlETLegendmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TLegend*>*)
   {
      vector<TLegend*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TLegend*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TLegend*>", -2, "vector", 469,
                  typeid(vector<TLegend*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETLegendmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TLegend*>) );
      instance.SetNew(&new_vectorlETLegendmUgR);
      instance.SetNewArray(&newArray_vectorlETLegendmUgR);
      instance.SetDelete(&delete_vectorlETLegendmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETLegendmUgR);
      instance.SetDestructor(&destruct_vectorlETLegendmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TLegend*> >()));

      ::ROOT::AddClassAlternate("vector<TLegend*>","std::__1::vector<TLegend*, std::__1::allocator<TLegend*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TLegend*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETLegendmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TLegend*>*)0x0)->GetClass();
      vectorlETLegendmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETLegendmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETLegendmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLegend*> : new vector<TLegend*>;
   }
   static void *newArray_vectorlETLegendmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLegend*>[nElements] : new vector<TLegend*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETLegendmUgR(void *p) {
      delete ((vector<TLegend*>*)p);
   }
   static void deleteArray_vectorlETLegendmUgR(void *p) {
      delete [] ((vector<TLegend*>*)p);
   }
   static void destruct_vectorlETLegendmUgR(void *p) {
      typedef vector<TLegend*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TLegend*>

namespace ROOT {
   static TClass *vectorlETH1mUgR_Dictionary();
   static void vectorlETH1mUgR_TClassManip(TClass*);
   static void *new_vectorlETH1mUgR(void *p = 0);
   static void *newArray_vectorlETH1mUgR(Long_t size, void *p);
   static void delete_vectorlETH1mUgR(void *p);
   static void deleteArray_vectorlETH1mUgR(void *p);
   static void destruct_vectorlETH1mUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH1*>*)
   {
      vector<TH1*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH1*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH1*>", -2, "vector", 469,
                  typeid(vector<TH1*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH1mUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TH1*>) );
      instance.SetNew(&new_vectorlETH1mUgR);
      instance.SetNewArray(&newArray_vectorlETH1mUgR);
      instance.SetDelete(&delete_vectorlETH1mUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH1mUgR);
      instance.SetDestructor(&destruct_vectorlETH1mUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH1*> >()));

      ::ROOT::AddClassAlternate("vector<TH1*>","std::__1::vector<TH1*, std::__1::allocator<TH1*> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TH1*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH1mUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TH1*>*)0x0)->GetClass();
      vectorlETH1mUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH1mUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH1mUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH1*> : new vector<TH1*>;
   }
   static void *newArray_vectorlETH1mUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH1*>[nElements] : new vector<TH1*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH1mUgR(void *p) {
      delete ((vector<TH1*>*)p);
   }
   static void deleteArray_vectorlETH1mUgR(void *p) {
      delete [] ((vector<TH1*>*)p);
   }
   static void destruct_vectorlETH1mUgR(void *p) {
      typedef vector<TH1*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TH1*>

namespace ROOT {
   static TClass *maplEstringcOvectorlETLegendmUgRsPgR_Dictionary();
   static void maplEstringcOvectorlETLegendmUgRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlETLegendmUgRsPgR(void *p = 0);
   static void *newArray_maplEstringcOvectorlETLegendmUgRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlETLegendmUgRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlETLegendmUgRsPgR(void *p);
   static void destruct_maplEstringcOvectorlETLegendmUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<TLegend*> >*)
   {
      map<string,vector<TLegend*> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<TLegend*> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<TLegend*> >", -2, "map", 898,
                  typeid(map<string,vector<TLegend*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlETLegendmUgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,vector<TLegend*> >) );
      instance.SetNew(&new_maplEstringcOvectorlETLegendmUgRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlETLegendmUgRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlETLegendmUgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlETLegendmUgRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlETLegendmUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<TLegend*> > >()));

      ::ROOT::AddClassAlternate("map<string,vector<TLegend*> >","std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> >, std::__1::vector<TLegend*, std::__1::allocator<TLegend*> >, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > >, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > const, std::__1::vector<TLegend*, std::__1::allocator<TLegend*> > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,vector<TLegend*> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlETLegendmUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,vector<TLegend*> >*)0x0)->GetClass();
      maplEstringcOvectorlETLegendmUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlETLegendmUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlETLegendmUgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<TLegend*> > : new map<string,vector<TLegend*> >;
   }
   static void *newArray_maplEstringcOvectorlETLegendmUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<TLegend*> >[nElements] : new map<string,vector<TLegend*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlETLegendmUgRsPgR(void *p) {
      delete ((map<string,vector<TLegend*> >*)p);
   }
   static void deleteArray_maplEstringcOvectorlETLegendmUgRsPgR(void *p) {
      delete [] ((map<string,vector<TLegend*> >*)p);
   }
   static void destruct_maplEstringcOvectorlETLegendmUgRsPgR(void *p) {
      typedef map<string,vector<TLegend*> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<TLegend*> >

namespace ROOT {
   static TClass *maplEstringcOvectorlETH1mUgRsPgR_Dictionary();
   static void maplEstringcOvectorlETH1mUgRsPgR_TClassManip(TClass*);
   static void *new_maplEstringcOvectorlETH1mUgRsPgR(void *p = 0);
   static void *newArray_maplEstringcOvectorlETH1mUgRsPgR(Long_t size, void *p);
   static void delete_maplEstringcOvectorlETH1mUgRsPgR(void *p);
   static void deleteArray_maplEstringcOvectorlETH1mUgRsPgR(void *p);
   static void destruct_maplEstringcOvectorlETH1mUgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<string,vector<TH1*> >*)
   {
      map<string,vector<TH1*> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<string,vector<TH1*> >));
      static ::ROOT::TGenericClassInfo 
         instance("map<string,vector<TH1*> >", -2, "map", 898,
                  typeid(map<string,vector<TH1*> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEstringcOvectorlETH1mUgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(map<string,vector<TH1*> >) );
      instance.SetNew(&new_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetNewArray(&newArray_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetDelete(&delete_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetDeleteArray(&deleteArray_maplEstringcOvectorlETH1mUgRsPgR);
      instance.SetDestructor(&destruct_maplEstringcOvectorlETH1mUgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<string,vector<TH1*> > >()));

      ::ROOT::AddClassAlternate("map<string,vector<TH1*> >","std::__1::map<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> >, std::__1::vector<TH1*, std::__1::allocator<TH1*> >, std::__1::less<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > >, std::__1::allocator<std::__1::pair<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > const, std::__1::vector<TH1*, std::__1::allocator<TH1*> > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<string,vector<TH1*> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEstringcOvectorlETH1mUgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<string,vector<TH1*> >*)0x0)->GetClass();
      maplEstringcOvectorlETH1mUgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEstringcOvectorlETH1mUgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<TH1*> > : new map<string,vector<TH1*> >;
   }
   static void *newArray_maplEstringcOvectorlETH1mUgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<string,vector<TH1*> >[nElements] : new map<string,vector<TH1*> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      delete ((map<string,vector<TH1*> >*)p);
   }
   static void deleteArray_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      delete [] ((map<string,vector<TH1*> >*)p);
   }
   static void destruct_maplEstringcOvectorlETH1mUgRsPgR(void *p) {
      typedef map<string,vector<TH1*> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<string,vector<TH1*> >

namespace {
  void TriggerDictionaryInitialization_MHists_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./MHists.C",
0
    };
    static const char* includePaths[] = {
"/Users/arkasantra/Root6_Install/include",
"/Users/arkasantra/Root6_Install/etc/",
"/Users/arkasantra/Root6_Install/etc//cling",
"/Users/arkasantra/Root6_Install/include/",
"/Users/arkasantra/Root6_Install/include",
"/Users/arkasantra/Root6_Install/include/",
"/Users/arkasantra/arka/Sasha_Work/lxmacro/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MHists_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(file_name@@@/Users/arkasantra/arka/Sasha_Work/lxmacro/./MHists.h)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$MHists.h")))  __attribute__((annotate("$clingAutoload$./MHists.C")))  MHists;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MHists_C_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./MHists.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"MHist_test_v1", payloadCode, "@",
"MHists", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MHists_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MHists_C_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MHists_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MHists_C_ACLiC_dict() {
  TriggerDictionaryInitialization_MHists_C_ACLiC_dict_Impl();
}
