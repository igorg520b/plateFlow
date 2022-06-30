#include "meshrepresentation.h"
#include "mesh.h"
#include "element.h"
#include "cohesivezone.h"
#include "model.h"

#include "spdlog/spdlog.h"

icy::MeshRepresentation::MeshRepresentation()
{
    int nLut = sizeof lutArrayTemperatureAdj / sizeof lutArrayTemperatureAdj[0];
    hueLut->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        hueLut->SetTableValue(i, lutArrayTemperatureAdj[i][0],
                              lutArrayTemperatureAdj[i][1],
                              lutArrayTemperatureAdj[i][2], 1.0);

    nLut = sizeof lutArrayPastel / sizeof lutArrayPastel[0];
    hueLut_pastel->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
        hueLut_pastel->SetTableValue(i, lutArrayPastel[i][0],
                              lutArrayPastel[i][1],
                              lutArrayPastel[i][2], 1.0);

    hueLut_czs->SetNumberOfTableValues(3);
    hueLut_czs->SetTableValue(0, 0.95, 0.97, 0.91, 0.05);
    hueLut_czs->SetTableValue(1, 0.25, 0.87, 0.21, 0.5);
    hueLut_czs->SetTableValue(2, 0.95, 0.37, 0.31, 0.9);

    // initialize various VTK objects
    visualized_values->SetName("visualized_values");
    visualized_values_czs->SetName("visualized_values_czs");

    ugrid_deformable->SetPoints(points_deformable);
    dataSetMapper_deformable->SetInputData(ugrid_deformable);
    dataSetMapper_deformable->UseLookupTableScalarRangeOn();
//    dataSetMapper_deformable->SetLookupTable(hueLut);

    actor_mesh_deformable->SetMapper(dataSetMapper_deformable);
    actor_mesh_deformable->GetProperty()->VertexVisibilityOff();
    actor_mesh_deformable->GetProperty()->EdgeVisibilityOn();
    actor_mesh_deformable->GetProperty()->SetColor(0.84938, 0.872213, 0.848103);
    actor_mesh_deformable->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
    actor_mesh_deformable->PickableOff();
    actor_mesh_deformable->GetProperty()->SetLineWidth(0.6);
    actor_mesh_deformable->GetProperty()->SetRenderLinesAsTubes(true);

//    actor_mesh_deformable->GetProperty()->LightingOff();
//    actor_mesh_deformable->GetProperty()->ShadingOff();
//    actor_mesh_deformable->GetProperty()->SetInterpolationToFlat();

    // indenter
    cylinder->SetResolution(50);
    cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
    cylinderActor->SetMapper(cylinderMapper);
    cylinderActor->GetProperty()->SetColor(0.14,0.15,0.4);
    cylinderActor->GetProperty()->SetRepresentationToWireframe();
//    cylinderActor->RotateX(30.0);
//    cylinderActor->RotateY(-45.0);

    // CZs
    ugrid_czs->SetPoints(points_deformable);
    mapper_czs->SetInputData(ugrid_czs);
    actor_czs->PickableOff();
    actor_czs->SetMapper(mapper_czs);
    actor_czs->GetProperty()->EdgeVisibilityOff();
//    actor_czs->GetProperty()->SetColor(0.54938, 0.572213, 0.548103);
//    actor_czs->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
//    actor_czs->GetProperty()->SetLineWidth(0.3);

    ugrid_czs->GetCellData()->AddArray(visualized_values_czs);
    ugrid_czs->GetCellData()->SetActiveScalars("visualized_values_czs");

    mapper_czs->UseLookupTableScalarRangeOn();
    mapper_czs->SetLookupTable(hueLut_czs);
    mapper_czs->SetScalarModeToUseCellData();
    mapper_czs->ScalarVisibilityOn();
}

void icy::MeshRepresentation::UnsafeSynchronizePositions()
{

    mesh->mutexMeshUpdate.lock();
    for(icy::Node *nd : mesh->nodes)
        points_deformable->SetPoint((vtkIdType)nd->globId, nd->xn.data());
    points_deformable->Modified();
    mesh->mutexMeshUpdate.unlock();

    double indenterCenter = 1. + model->prms.p.R_indenter - model->indenter_position;
    cylinder->SetCenter(0.0, 0.0, indenterCenter);
    cylinder->SetRadius(model->prms.p.R_indenter);
    cylinder->SetHeight(1.0);
}


void icy::MeshRepresentation::UnsafeSynchronizeTopology()
{
    mesh->mutexMeshUpdate.lock();

    points_deformable->SetNumberOfPoints(mesh->nodes.size());
    for(icy::Node *nd : mesh->nodes)
    {
        double x[3] {nd->xn[0],nd->xn[1],nd->xn[2]};
        points_deformable->SetPoint((vtkIdType)nd->globId, x);
    }

    // deformable material - elements
    cellArray_deformable->Reset();
    for(icy::Element *tr : mesh->elems)
    {
        vtkIdType pts[4];
        for(int k=0;k<4;k++) pts[k] = tr->nds[k]->globId;
        cellArray_deformable->InsertNextCell(4, pts);
    }
    ugrid_deformable->SetCells(VTK_TETRA, cellArray_deformable);
//    ugrid_deformable->Modified();

    // CZs
    cellArray_czs->Reset();
    for(icy::CohesiveZone *cz : mesh->czs)
    {
        vtkIdType pts[3] {cz->nds[0]->globId,cz->nds[1]->globId,cz->nds[2]->globId};
        cellArray_czs->InsertNextCell(3,pts);
    }
    ugrid_czs->SetCells(VTK_TRIANGLE, cellArray_czs);
    visualized_values_czs->SetNumberOfValues(mesh->czs.size());

    mesh->mutexMeshUpdate.unlock();
    UnsafeSynchronizeValues();
}


void icy::MeshRepresentation::UnsafeSynchronizeValues()
{
    if(mesh->nodes.size()==0)
    {
        visualized_values->SetNumberOfValues(0);
        return;
    }

    mesh->mutexMeshUpdate.lock();

    int nLut = sizeof lutArrayPastel / sizeof lutArrayPastel[0];
    double minmax[2];
    switch(VisualizingVariable)
    {
    case VisOpt::none:
        break;
    case VisOpt::grain_id:
        visualized_values->SetNumberOfValues(mesh->elems.size());
        for(std::size_t i=0;i<mesh->elems.size();i++) visualized_values->SetValue(i, (double)(mesh->elems[i]->grainId%nLut));
        hueLut_pastel->SetTableRange(0, nLut-1);
        break;

    case VisOpt::stress_hydrostatic:
        visualized_values->SetNumberOfValues(mesh->elems.size());
        for(std::size_t i=0;i<mesh->elems.size();i++)
            visualized_values->SetValue(i, mesh->elems[i]->hydrostatic_stress);
        visualized_values->Modified();
        visualized_values->GetValueRange(minmax);
        hueLut->SetTableRange(minmax[0], minmax[1]);
        break;

    case VisOpt::stress_p1:
        visualized_values->SetNumberOfValues(mesh->elems.size());
        for(std::size_t i=0;i<mesh->elems.size();i++)
            visualized_values->SetValue(i, mesh->elems[i]->principal_stress[0]);
        visualized_values->Modified();
        visualized_values->GetValueRange(minmax);
        hueLut->SetTableRange(minmax[0], minmax[1]);
        break;

    case VisOpt::energy_density:
        visualized_values->SetNumberOfValues(mesh->elems.size());
        for(std::size_t i=0;i<mesh->elems.size();i++)
            visualized_values->SetValue(i, mesh->elems[i]->strain_energy_density);
        visualized_values->Modified();
        visualized_values->GetValueRange(minmax);
        hueLut->SetTableRange(minmax[0], minmax[1]);
        break;


    default:
        break;
    }

    for(std::size_t i=0;i<mesh->czs.size();i++)
    {
        icy::CohesiveZone *cz = mesh->czs[i];
        visualized_values_czs->SetValue(i, cz->isActive ? (cz->isDamaged ? 1 : 0) : 2);
    }
    visualized_values_czs->Modified();
    hueLut_czs->SetTableRange(0,2);


    mesh->mutexMeshUpdate.unlock();
}


void icy::MeshRepresentation::ChangeVisualizationOption(int option)
{
    spdlog::info("icy::Model::ChangeVisualizationOption {}", option);
    VisualizingVariable = (VisOpt)option;

    ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
    ugrid_deformable->GetCellData()->RemoveArray("visualized_values");

    switch(VisualizingVariable)
    {
    case VisOpt::grain_id:
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->AddArray(visualized_values);
        ugrid_deformable->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUseCellData();
        dataSetMapper_deformable->ScalarVisibilityOn();
        dataSetMapper_deformable->SetLookupTable(hueLut_pastel);
        break;

    case VisOpt::energy_density:
    case VisOpt::stress_hydrostatic:
    case VisOpt::stress_p1:
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->AddArray(visualized_values);
        ugrid_deformable->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUseCellData();
        dataSetMapper_deformable->ScalarVisibilityOn();
        dataSetMapper_deformable->SetLookupTable(hueLut);
        break;



/*
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        ugrid_deformable->GetPointData()->AddArray(visualized_values);
        ugrid_deformable->GetPointData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUsePointData();
        dataSetMapper_deformable->ScalarVisibilityOn();
        break;
*/
    }
    UnsafeSynchronizeValues();

}
