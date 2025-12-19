/*
# Copyright (C) 2021 Sophie Bini, Gabriele Vedovato
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#define XIFO 4

#pragma GCC system_header

#include "cwb.hh"
#include <python3.10/Python.h>


// ---------------------------------------------------------------------------------
// WHAT IS?
// this plugin compute the gitchiness uging the autoencoder
// HOW TO CONFIGURE THE AE PLUGIN
// the following is an example : must be included in the config/user_parameters.C file
// see the 'AUTOENCODER DEFAULT USER DEFINES' function for the full description of parameters
// ---------------------------------------------------------------------------------
/*
  plugin = TMacro(gSystem->ExpandPathName("$HOME_CWB/plugins/CWB_Plugin_AE.C"));        // Macro source

  TString optae = "";                           			// NOTE : add space at the end of each line
  optae += "ae_weights_fname=/user/albert.einstein/ae_weights.h5 ";     // autoencoder weights file name
  optae += "ae_output_root=false ";                   			// disable output to root file by the ae plugin

  strcpy(parPlugin,optae.Data());               			// set AE plugin parameters
  strcpy(comment,"ae configuration example");
*/

// ---------------------------------------------------------------------------------
// AUTOENCODER DEFINES
// ---------------------------------------------------------------------------------

#define AE_DIM			1	// output root autoencoder array

// ---------------------------------------------------------------------------------
// AUTOENCODER DEFAULT USER DEFINES
// ---------------------------------------------------------------------------------

#define AE_ENABLE        true    		// if true ->   enable autoencoder
#define AE_PYTHON_PATH   "$HOME_CWB/scripts"    // autoencoder python path name
#define AE_WEIGHTS_FNAME ""      		// autoencoder weights file name
#define AE_OUTPUT_ROOT   false   		// (def=false) set to false if this plugin is used with QLveto (QLveto plugin output par to output root file)
                         			// NOTE: QLveto plugin must be declared as the last one in cwb_mplugin

// ---------------------------------------------------------------------------------
// USER CONFIG OPTIONS
// ---------------------------------------------------------------------------------

struct uoptions {
  bool    enable;
  TString python_path;
  TString weights_fname;
  bool    output_root;
};

// ---------------------------------------------------------------------------------
// Embedded Python variables used to call autoencoder python script
// ---------------------------------------------------------------------------------

struct PY {
  PyObject *pName;
  PyObject *pModule;
  PyObject *pDict;
  PyObject *pClass;
  PyObject *pInstance;
  PyObject *pValue;
};

// ---------------------------------------------------------------------------------
// Global Variables
// ---------------------------------------------------------------------------------

uoptions   gOPT;	// global User Options
PY	   gPY;		// global autoencoder embedded python variables
TTree*     gTREE;       // output tree file name
TString    gOUTPUT;     // output root file name
float	   gAE[AE_DIM]; // autoencoder output values array

// ---------------------------------------------------------------------------------
// FUNCTIONS
// ---------------------------------------------------------------------------------

void  ResetUserOptions();
void  ReadUserOptions(TString options);
void  PrintUserOptions(CWB::config* cfg);

void  InitializeAutoEncoder();
float GetGlitchiness(wavearray<double>* w);
void  FinalizeAutoEncoder();

void  ClearWaveforms(detector* ifo);
void  SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_ae);
void  DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor);

void 
CWB_Plugin(TFile* jfile, CWB::config* cfg, network* NET, WSeries<double>* x, TString ifo, int type)  {
//!MISCELLANEA
// Extract whitened reconstructed waveforms, and compute the autoencoder output values

//  cout << endl;
//  cout << "-----> CWB_Plugin_AE.C" << endl;
//  cout << "ifo " << ifo.Data() << endl;
//  cout << "type " << type << endl;
//  cout << endl;

  if(type==CWB_PLUGIN_CONFIG) {  
    cout << endl;
    cout << "-----> CWB_Plugin_AE.C" << endl;
    cout << endl;

    ResetUserOptions();                                 // set default config options
    ReadUserOptions(cfg->parPlugin);                    // user config options : read from parPlugin

    if(gOPT.enable) {
      PrintUserOptions(cfg);               		// print config options
      InitializeAutoEncoder();     			// init python autoencoder class
      CWB::Toolbox::checkFile(gOPT.weights_fname);	// check if autoencoder weights file exists
      cfg->outPlugin=true;  				// disable built-in output root file
    }
  }

  if(type==CWB_PLUGIN_ILIKELIHOOD) {                    // DEFINE WAVE TREE - FIX CASE WHEN NO EVENTS HAVE BEEN FOUND !!!
    if(!gOPT.enable) return;
    netevent* EVT;
    SetOutputFile(NET, EVT, cfg, false);                // set default output waveburst tree in root file
  }

  if(type==CWB_PLUGIN_OLIKELIHOOD) {

    if(!gOPT.enable) return;

    // import ifactor
    int gIFACTOR=-1; IMPORT(int,gIFACTOR)
    double factor = cfg->simulation==3||cfg->simulation==4 ? -gIFACTOR : cfg->factors[gIFACTOR];
    double ofactor=0;
    if(cfg->simulation==4)      ofactor=-gIFACTOR;
    else if(cfg->simulation==3) ofactor=-gIFACTOR;
    else                        ofactor=factor;

    int nIFO = NET->ifoListSize();                      // number of detectors
    int K = NET->nLag;                                  // number of time lag
    int rate = 0;                                       // select all resolutions
    netevent* EVT;
    wavearray<double> id;

    SetOutputFile(NET, EVT, cfg, false);                // set output root file

    for(int k=0; k<K; k++) {                            // loop over the lags

      id = NET->getwc(k)->get(const_cast<char*>("ID"), 0, 'L', rate);

      for(int j=0; j<(int)id.size(); j++) {             // loop over cluster index

        int ID = size_t(id.data[j]+0.5);

        if(NET->getwc(k)->sCuts[ID-1]!=-1)  continue;   // skip rejected/processed clusters

        EVT->output2G(NULL,NET,ID,k,ofactor);           // get reconstructed parameters

        // extract whitened reconstructed waveforms
        // compute ae[0] using the weighted of glness values
        gAE[0]=0.;
        for(int n=0; n<nIFO; n++) {
	   // reconstructed whitened waveform
           NET->getMRAwave(ID,k,'S',0,true);
           detector* pD = NET->getifo(n);
	   gAE[0]+=EVT->sSNR[n]*GetGlitchiness(&(pD->waveForm));
        }
        gAE[0]/=EVT->likelihood;
        cout << endl << "autoencoder   :" << " ae[0] = " << gAE[0] << endl << endl;

        SetOutputFile(NET, EVT, cfg, true);             // set output root file

        if(gOPT.output_root) {
          DumpOutputFile(NET, EVT, cfg, ID, k, ofactor);        // dump event to output root file

          for(int n=0;n<nIFO;n++) {
            detector* pD = NET->getifo(n);
            if(!cfg->simulation) ClearWaveforms(pD);            // release waveform memory
          }
        }
      }
    }

    jfile->cd();
    if(EVT) delete EVT;
  }
  return;
}

void SetOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, bool dump_ae) {

   // import slagShift
   float* gSLAGSHIFT=NULL; IMPORT(float*,gSLAGSHIFT)

   int nIFO = NET->ifoListSize();                       // number of detectors

   // search output root file in the system list
   TFile* froot = NULL;
   TList *files = (TList*)gROOT->GetListOfFiles();
   gOUTPUT="";
   if (files) {
     TIter next(files);
     TSystemFile *file;
     TString fname;
     bool check=false;
     while ((file=(TSystemFile*)next())) {
        fname = file->GetName();
        // set output root file as the current file
        if(fname.Contains("wave_")) {
          froot=(TFile*)file;froot->cd();
          gOUTPUT=fname;
          gOUTPUT.ReplaceAll(".root.tmp",".txt");
          //cout << "output file name : " << fname << endl;
        }
     }
     if(!froot) {
       cout << "CWB_Plugin_AE.C : Error - output root file not found" << endl;
       gSystem->Exit(1);
     }
   } else {
     cout << "CWB_Plugin_AE.C : Error - output root file not found" << endl;
     gSystem->Exit(1);
   }

   gTREE = (TTree *) froot->Get("waveburst");
   if(gTREE!=NULL) {
     EVT = new netevent(gTREE,nIFO);
     if(dump_ae) {
       TBranch* branch;
       bool ae_exists=false;
       TIter next(gTREE->GetListOfBranches());
       while ((branch=(TBranch*)next())) {
         if(TString("ae").CompareTo(branch->GetName())==0) ae_exists=true;
       }
       next.Reset();
       if(!ae_exists) gTREE->Branch("ae", gAE, TString::Format("ae[%i]/F",AE_DIM));
     }
   } else {
     EVT = new netevent(nIFO);
     gTREE = EVT->setTree();
   }
   EVT->setSLags(gSLAGSHIFT);                          // set slags into netevent           
   EVT->Psave=cfg->Psave;
}

void DumpOutputFile(network* NET, netevent* &EVT, CWB::config* cfg, int ID, int k, int factor) {

   if(cfg->dump) EVT->dopen(gOUTPUT.Data(),const_cast<char*>("a"),false);
   EVT->output2G(gTREE,NET,ID,k,factor);             // get reconstructed parameters
   if(cfg->dump) EVT->dclose();
   if(!cfg->cedDump) NET->getwc(k)->sCuts[ID-1]=1;   // mark as processed
}

void
ClearWaveforms(detector* ifo) {

  int n;

  n = ifo->IWFP.size();
  for (int i=0;i<n;i++) {
    wavearray<double>* wf = (wavearray<double>*)ifo->IWFP[i];
    delete wf;
  }
  ifo->IWFP.clear();
  ifo->IWFID.clear();

  n = ifo->RWFP.size();
  for (int i=0;i<n;i++) {
    wavearray<double>* wf = (wavearray<double>*)ifo->RWFP[i];
    delete wf;
  }
  ifo->RWFP.clear();
  ifo->RWFID.clear();
}

void
InitializeAutoEncoder() {

  gSystem->Setenv("TF_CPP_MIN_LOG_LEVEL","2");	// INFO and WARNING and ERROR tensorflow messages are not printed

  gSystem->Setenv("PYTHONDONTWRITEBYTECODE","x cwb_autoencoder"); // disable __pycache__/*.pyc

  Py_Initialize();

  // initialize dummy argc,argv
  wchar_t const *dummy_args[] = {L"Python", NULL};  // const is needed because literals must not be modified
  wchar_t const **argv = dummy_args;
  int             argc = sizeof(dummy_args)/sizeof(dummy_args[0])-1;
  PySys_SetArgv(argc, const_cast<wchar_t **>(argv)); // const_cast allowed, because PySys_SetArgv doesn't change arg

  // set module path
  TString module_path = gSystem->ExpandPathName(gOPT.python_path.Data());
  PyObject *sys_path = PySys_GetObject("path");
  PyList_Append(sys_path, PyUnicode_FromString(module_path.Data()));

  gPY.pName = PyUnicode_FromString((char*)"cwb_autoencoder");
  gPY.pModule = PyImport_Import(gPY.pName);
  Py_DECREF(gPY.pName);
  if(gPY.pModule == NULL) {PyErr_Print();exit(1);}

  gPY.pDict = PyModule_GetDict(gPY.pModule);
  if(gPY.pDict == NULL) {PyErr_Print();exit(1);}

  // Build the name of a callable class 
  gPY.pClass = PyDict_GetItemString(gPY.pDict, "AutoEncoder");
  if(gPY.pClass == NULL) {PyErr_Print();exit(1);}

  // Create an instance of the class
  if(PyCallable_Check(gPY.pClass)) gPY.pInstance = PyObject_CallObject(gPY.pClass, NULL);

  // set autoencoder weights
  PyObject* weights_fname = PyUnicode_FromString(gOPT.weights_fname);
  PyObject* pValue = PyObject_CallMethod(gPY.pInstance, (char*)"set_weights", (char*)"(O)", weights_fname);

  return;
}

float
GetGlitchiness(wavearray<double>* w) {

  float rate = w->rate();
  // w data are copied to string
  std::string sdata = "";
  for(int i=0;i<w->size();i++) sdata+=TString::Format(" %g",w->data[i]);
  // call python method glness
  PyObject* data = PyUnicode_FromString(sdata.c_str());
  PyObject* pValue = PyObject_CallMethod(gPY.pInstance, (char*)"get_glness", (char*)"(Of)", data, rate);

  // read glness
  float glness=-1;
  if(pValue != NULL) {
     glness=PyFloat_AsDouble(pValue);
     //printf("Return of call : %g\n", glness);
     Py_DECREF(pValue);
  } else {
     PyErr_Print();
     exit(1);
  }

  return glness;
}

void
FinalizeAutoEncoder() {

  // Clean up
  Py_DECREF(gPY.pModule);
  Py_DECREF(gPY.pName);
  return;
}

void ResetUserOptions() {

  gOPT.enable    	= AE_ENABLE;
  gOPT.python_path   	= AE_PYTHON_PATH;
  gOPT.weights_fname   	= AE_WEIGHTS_FNAME;
  gOPT.output_root  	= AE_OUTPUT_ROOT;
}

void ReadUserOptions(TString options) {

  if(options.CompareTo("")!=0) {
    cout << options << endl;
    if(!options.Contains("--")) {  // parameters are used only by cwb_inet

      TObjArray* token = TString(options).Tokenize(TString(' '));
        for(int j=0;j<token->GetEntries();j++){

        TObjString* tok = (TObjString*)token->At(j);
        TString stok = tok->GetString();

        if(stok.Contains("ae_enable=")) {
          TString enable=stok;
          enable.Remove(0,enable.Last('=')+1);
          gOPT.enable=(enable=="true")?true:false;
        }

        if(stok.Contains("ae_python_path=")) {
          TString python_path=stok;
          python_path.Remove(0,python_path.Last('=')+1);
          gOPT.python_path=python_path;
        }

        if(stok.Contains("ae_weights_fname=")) {
          TString weights_fname=stok;
          weights_fname.Remove(0,weights_fname.Last('=')+1);
          gOPT.weights_fname=weights_fname;
        }

        if(stok.Contains("ae_output_root=")) {
          TString output_root=stok;
          output_root.Remove(0,output_root.Last('=')+1);
          gOPT.output_root=(output_root=="true")?true:false;
        }
      }
    }
  }
}

void PrintUserOptions(CWB::config* cfg) {

  cout << "-----------------------------------------"     << endl;
  cout << "AE config options                        "     << endl;
  cout << "-----------------------------------------"     << endl << endl;
  cout << "AE_ENABLE              " << gOPT.enable        << endl;
  cout << "AE_PYTHON_PATH         " << gOPT.python_path   << endl;
  cout << "AE_WEIGHTS_FNAME       " << gOPT.weights_fname << endl;
  cout << "AE_OUTPUT_ROOT         " << gOPT.output_root   << endl;

  cout << endl;
}

