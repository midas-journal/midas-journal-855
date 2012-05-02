//
//  DataWrappers.h
//  
//
//  Created by Mathieu De Craene on 2/4/12.
//  Copyright 2012 UPF. All rights reserved.
//

#ifndef _DataWrappers_h
#define _DataWrappers_h

#include "itkProcessObject.h"
#include "itkDataObject.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"


class ITK_EXPORT MeshContainer : public itk::DataObject
{
public:
	//! Standard itk typedefs
	typedef MeshContainer Self;
	typedef itk::DataObject ParentType;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self>  ConstPointer;
    
	/** Run-time type information (and related methods). */
	itkTypeMacro(Self,ParentType);
	itkFactorylessNewMacro(Self);
    
	void AddInput( vtkPolyData* pd )
	{
        pds.push_back(pd);
	}
    
    
	//! Holds the meshes
	std::vector<vtkSmartPointer<vtkPolyData> > pds;
};

#endif
