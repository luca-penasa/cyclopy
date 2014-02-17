# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np
import progressbar as pbar
from scipy.optimize import leastsq

pb = pbar.ProgressBar()


class UndeformedModel:
    '''
    Given a point_list of type [llUnstructuredGrid_1, llUnstructuredGrid_2, ...]
    A list is used couse in some cases we are working with cloud that were registered but not merged.
    In each of these clouds there are points with x,y,z,class_id. class_id is an array of the vtkUnstructuredGrid object.
    the class id array provide a subdivision of points into different classes
    Every class should have a unique W value, that is unknown: point into a class represent a stratigraphic bedding plane,
    so the stratigraphic position (SP) is unique for that class. 
    this class computes parameters for function of type W(x,y,z) = a*x + b*y + c*z , 
    and estimates the SP for each class (bedding plane)
    with the constrain that |grad(W)|_2 = 1. Then the problem is non-linear. Levenberg-Marquard algorithm is used.
    Once the model is parametrized We can ask W for a point simply calling the class UndeformedModel(point)
    and function gradient (aka stratigraphic normal - that is constant - using GetWGrad(point) function)
    
    Note: similar results could probably be obtained using OLS (orthogonal least squares) but this method would always need 
    to know the W value for each class: this situation exists only if we would use previously measured stratigrapohic sections
    so to know the strata thickness.
    
    Many functions are internal TODO: rewrite with lowercase!
    TO BE USED:
    - SetInput -> as input a numpy array as from pickedpointshandler
    - GetOutput -> a numpy array with points and residuals [x,y,z,residuals] <- vertical vectors
    '''
    
    def __init__(self):
        #for representing the plane
        self.normal = np.array([0,0,1])
        self.data_changed = False
        
    def EvaluateFunction(self, point):
        if len(point.shape) == 1:
            point = np.array([point])
        self.Update()
        a,b,c = self.normal
        x,y,z = point.T
        result = a*x + b*y + c*z
        return result
    
    def EvaluateGradient(self, point=None):
        self.Update()
        if point != None:
            shp = point.shape
            if len(shp) == 2:
                result = np.asarray([self.normal for i in range(shp[0])])
            else:
                result = self.normal
        else:
            result = self.normal
        return result
    
    def SetInput(self, points):
        self.points = points[:,0:3]
        self.class_id = points[:,3]
        self.data_changed = True
          
    def Update(self):
        if self.data_changed == True:
            self.N_points = len(self.points)
            self.N_classes = len(np.unique(self.class_id))
            self.unique_class_ids = np.unique(self.class_id).astype(int)
            self.UpdateWForClasses()
            self.Refine()
            self.data_changed = False
        
    def GetPointsInClass(self, class_index):
        points_to_take = self.class_id == class_index
        points = self.points[points_to_take]
        if len(points) == 0:
            return None
        return points
    
    def SetPlanePars(self, a,b,c):
        normal = np.array([a,b,c])
        #we have a constrain:
        #normal need to be normalized -> function have ||gradient||_2 = 1
        normal = normal / np.linalg.norm(normal)
        self.normal = normal
        self.UpdateWForClasses()
        
    def GetWForPointsInClass(self, class_id):
        a,b,c = self.normal    
        points = self.GetPointsInClass(class_id).T
        if len(points) == 1:
            x,y,z = points[0]
        else:
            x,y,z = points
        ws = a*x + b*y + c*z
        return ws
        
    def UpdateWForClasses(self):
        class_w = []
        for cid in self.unique_class_ids:
            this_w = np.mean(self.GetWForPointsInClass(cid))
            class_w.append(this_w)
            
        self.classes_w = np.array(class_w)
    
    def GetPredictionErrorForClass(self, class_id):
        points_w = self.GetWForPointsInClass(class_id)
#        class_w = self.classes_w[class_id]
        class_w = self.GetWForClass(class_id)
        error = np.sqrt(np.sum((points_w - class_w)**2))
        return error
    
    def GetWForClass(self, cid):
        local_index = np.arange(self.N_classes)[self.unique_class_ids == cid]
        return self.classes_w[local_index]
    
    def GetPredictionErrorForAllClasses(self):
        errors = []
        for cid in self.unique_class_ids:
            class_error = self.GetPredictionErrorForClass(cid)
            errors.append(class_error)

        return np.array(errors)
    
    def GetResiduals(self, normal):
        a,b,c = normal
        self.SetPlanePars(a,b,c)
        residuals = self.GetPredictionErrorForAllClasses()
        return residuals
        
    def Refine(self):
        self.leastsq_solution = leastsq(self.GetResiduals, self.normal, ftol = 1e-6, full_output=False)
        
    def GetWForPoints(self):
        result = np.ones(len(self.GetPointsInClass(self.class_id[0])))
        result *= self.classes_w[0]
        for i,cid in enumerate(self.unique_class_ids[1::]):
            part = np.ones(len(self.GetPointsInClass(cid)))
            part *= self.classes_w[i]
            result = np.hstack([result, part])
        return result
        
#    def GetWForClass(self, class_id):
#        return self.classes_w[class_id]
    
    def GetOutput(self):
        self.Update() #automatically call the update method
        real_w = self.GetWForPoints()
        modeled_w = self.EvaluateFunction(self.points)
        output = real_w - modeled_w
        output = np.vstack([self.points.T, output]).T
        return output
        
    def GetComputedW(self):
        self.Update()
        self.Refine()
        return self.unique_class_ids, self.classes_w
    
    
class ResidualsModel:
    '''
    This class is used to achieve a perfect interpolation between picked points and the model. 
    Is based upon Radial Basis Function (RBF)
    Internal Functions are in lowercase
    '''
    def __init__(self):
        self.data_changed = False
        self.R = 500 #default radius
        self.smooth = 0
        
    def SetInput(self, points):
        self.x, self.y, self.z = points[:,0:3].T
        self.values = points[:,3]
        self.N = len(self.x)
        self.data_changed = True
    
    def SetRadius(self, R):
        self.R = R
        self.data_changed = True
        
    def SetSmooth(self, smooth):
        '''
        Remember that the smooth parameter have differents behavuiour depending on 
        the rbf used
        '''
        self.smooth = smooth
        self.data_changed = True
    
    def Update(self):
        if self.data_changed == True:
            self.A = self.compute_distance_matrix()
            self.nodes = self.compute_nodes()
            self.data_changed = False
        
    def compute_distance_matrix(self):
        A = np.zeros((self.N,self.N))
        for i in np.arange(self.N):
            p0 = np.array([self.x[i], self.y[i], self.z[i]])
            for j in np.arange(self.N):
                if i > j:
                    p = np.array([self.x[j], self.y[j], self.z[j]])
                    r = np.linalg.norm(p-p0)
                    A[i,j] = r
                    A[j,i] = r
                else:
                    pass
        return A
                    
    def thin_plate(self, r):
        result = r**2 * np.log(r)
        result[r == 0] = 0 # the spline is zero at zero
        return result
    
    
    def inverse_quadric(self,r):
        result = 1/(1 + (self.R*r)**2)
        return result
        
    def inverse_quadric_grad(self, delta, r):
        result = -2*delta*self.R**2 * 1 / (r**2*self.R**2 + 1)**2
        return result
    
    def elastic_def(self, r):
        result = (r**2/self.R**2) * np.log(r**2/self.R**2) + 1 - (r**2/self.R**2)
        result[r==0] = 1 # the limit of the function for r->0 is 1
        return result
    
    def elastic_def_grad(self, delta, r):
        result = 2*delta/self.R**2*np.log(r**2/self.R**2)
        result[r==0] = 0 #the function's limit for r->0 is 0
        result[delta==0] = 0 #as above
        return result
        
    def compute_nodes(self):
        '''
        computes nodes coefficients
        '''
        self.C = np.zeros((self.N, self.N))
        
        for i in np.arange(self.N):
            self.C[i] = self.elastic_def(self.A[i])
        if self.smooth != 0:
            self.C = self.C - np.eye(self.N) * self.smooth
        
        return  np.linalg.solve(self.C, self.values)
    
    
    def EvaluateFunction(self, points):
        self.Update()
        if len(points.shape) == 1:
            N_points = 1
            points = np.array([points])
        else:
            N_points = points.shape[0]
            
        nodes_xyz = np.vstack([self.x, self.y, self.z]).T
        
        
        distances = np.zeros((self.N, N_points))
        for i in np.arange(self.N):
            pn = nodes_xyz[i]
            for j in np.arange(N_points):
                p = points[j]
                distances[i,j] = np.sqrt(np.sum((p-pn)**2))
        
        
        rbf_value = self.elastic_def(distances)
        result = np.sum(self.nodes * rbf_value.T, axis = 1) 
        if len(result) == 1:
            result = result[0]
        return result
        
    def EvaluateGradient(self, points):
        self.Update()
        if len(points.shape) == 1:
            N_points = 1
            points = np.array([points])
        else:
            N_points = points.shape[0]
        
        
        nodes_xyz = np.vstack([self.x, self.y, self.z]).T #nodes positions
        
        #computing distances
        distances = np.zeros((self.N, N_points))
        for i in np.arange(self.N):
            pn = nodes_xyz[i]
            for j in np.arange(N_points):
                p = points[j]
                distances[i,j] = np.sqrt(np.sum((p-pn)**2))
        
        
        part_distances = np.zeros((self.N, N_points, 3))
        for i in np.arange(self.N):
            pn = nodes_xyz[i]
            for j in np.arange(N_points):
                p = points[j]
                part_distances[i,j] = p-pn
        
        
        rbf_value_x = self.elastic_def_grad(part_distances[:,:,0], distances)
        rbf_value_y = self.elastic_def_grad(part_distances[:,:,1], distances)
        rbf_value_z = self.elastic_def_grad(part_distances[:,:,2], distances)
        
        result_x = np.sum(self.nodes * rbf_value_x.T, axis = 1) 
        result_y = np.sum(self.nodes * rbf_value_y.T, axis = 1) 
        result_z = np.sum(self.nodes * rbf_value_z.T, axis = 1) 
        
        
        if len(result_x) == 1:
            result = np.array([result_x[0], result_y[0], result_z[0]]).T
        else:
            result = np.array([result_x, result_y, result_z]).T
    
        return result
        
        
class OutcropModel:
    def SetInput(self, points):
        self.umodel = UndeformedModel()
        self.umodel.SetInput(points)        
        self.rbfs = ResidualsModel()
        self.rbfs.SetInput(self.umodel.GetOutput())
        
        
    def EvaluateFunction(self, points, separated=True):
        umodel_w = self.umodel.EvaluateFunction(points)
        rbfs_w = self.rbfs.EvaluateFunction(points)
        if separated == True:
            return umodel_w, rbfs_w
        else:
            return umodel_w + rbfs_w
            
    def EvaluateFunction2(self, points, separated=True):
        N = len(points)
        umodel = []
        rbfs = []
        for point in pb(points):
            umodel_w = self.umodel.EvaluateFunction(point)
            rbfs_w = self.rbfs.EvaluateFunction(point)
            umodel.append(umodel_w)
            rbfs.append(rbfs_w)
        umodel = np.array(umodel)
        rbfs = np.array(rbfs)
            
        if separated == True:
            return umodel, rbfs
        else:
            return umodel + rbfs
        
    def EvaluateGradient(self, points):
        umodel_grad_w = self.umodel.EvaluateGradient(points)
        rbfs_grad_w = self.rbfs.EvaluateGradient(points)
        return umodel_grad_w + rbfs_grad_w
