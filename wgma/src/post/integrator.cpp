#include "post/integrator.hpp"

#include <TPZMaterial.h>
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <pzinterpolationspace.h>
#include <pzmultiphysicselement.h>
#include <pzvec_extras.h>

namespace wgma::post{

  Integrator::Integrator(TPZAutoPointer<TPZCompMesh> mesh,
                         std::set<int> matids,
                         int nThreads)
    : m_cmesh(mesh), m_nthreads(nThreads)
  {
    const auto meshdim = mesh->Dimension();
    
    if(matids.size() == 0){
      //insert every material of the mesh
      for(auto [id,matp] : mesh->MaterialVec()){
        auto bnd =
          dynamic_cast<TPZBndCond *>(matp);
        //we skip boundary materials
        if(matp && !bnd){
          matids.insert(id);
        }
      }
    }
    //now we get a list of all valid elements
    for(auto el : mesh->ElementVec()){
      if(!el || ! el->Reference()){ continue;}
      //checks if material is valid
      if(!el->HasMaterial(matids)){continue;}
      //checks whether the element has been refined
      const auto gel = el->Reference();
      if(gel->HasSubElement()){continue;}
      //checks element dimension
      const auto geldim = gel->Dimension();
      if(geldim != meshdim){continue;}
      AppendToVec(m_elvec,el);
    }
  }
  
  void Integrator::Integrate(const TPZVec<TPZCompEl*> &elvec)
  {
    auto mesh = this->Mesh();
    const auto meshdim = mesh->Dimension();
    

    //this lambda expression will be given to each thread
    auto ThreadWork = [this,&elvec](const int firstel,
                                    const int lastel,
                                    const int index){


      for(int iel = firstel; iel < lastel; iel++){
        auto el = elvec[iel];
        ElData data;
        InitData(el,data);
        auto &intrule = el->GetIntegrationRule();
        const int npts = intrule.NPoints();
        TPZManVector<REAL,3> pos(el->Dimension(),0);
        REAL weight;
        for(auto ipt = 0; ipt < npts; ipt++){
          intrule.Point(ipt, pos, weight);
          IntPointData(el,data,pos);
          Compute(data, weight, index);
        }
      }
    };

    //get the same number of threads from the assembly
    const int nthreads = this->NThreads();
    const int nel = elvec.size();
    if(nthreads < 1){
      ThreadWork(0,nel,0);
    }else{
      std::vector<std::thread> allthreads;
      //number of elements per thread (last one may have more elements)
      const int threadnel = nel / nthreads;
      for(int i = 0; i < nthreads; i++){
        int firstel = threadnel*i;
        int lastel = i == nthreads -1 ? nel : firstel + threadnel;
        allthreads.push_back(std::thread(ThreadWork,firstel,lastel,i));
      }
      for(int i = 0; i < nthreads; i++){
        allthreads[i].join();
      }
    }
  }


  
  void SingleSpaceIntegrator::InitData(TPZCompEl* el, ElData& data)
  {
    TPZMaterialDataT<CSTATE> &eldata = data;
    auto intel = dynamic_cast<TPZInterpolationSpace*>(el);
    intel->InitMaterialData(eldata);
  }
    
  void SingleSpaceIntegrator::IntPointData(TPZCompEl* el, ElData& data, TPZVec<REAL> &qsi)
  {
    TPZMaterialDataT<CSTATE> &eldata = data;
    eldata.fNeedsSol = true;
    auto intel = dynamic_cast<TPZInterpolationSpace*>(el);
		intel->ComputeRequiredData(eldata, qsi);
  }

  void MultiphysicsIntegrator::InitData(TPZCompEl* el, ElData& data)
  {
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(el);
    const int64_t nref = mfcel->NMeshes();
    TPZVec<TPZMaterialDataT<CSTATE>>& datavec = data;
    datavec.resize(nref);
    mfcel->InitMaterialData(datavec);
  }
    
  void MultiphysicsIntegrator::IntPointData(TPZCompEl* el, ElData& data, TPZVec<REAL> &qsi)
  {
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(el);
    const int64_t nref = mfcel->NMeshes();
    TPZVec<TPZMaterialDataT<CSTATE>>& datavec = data;
    for(int ir = 0; ir < nref; ir++){
      datavec[ir].fNeedsSol= true;
    }
    TPZManVector<TPZTransform<> > trvec;
    mfcel->AffineTransform(trvec);
    mfcel->ComputeRequiredData(qsi, trvec, datavec);
  }
  
};
