/**
 * \Generation of Spherocylindircal NP.
 */

//#define WARNINGS_ENABLED

#include <iostream>
#include <cstdlib>


#include "Include/tesselaSphere.h"

using namespace std;


struct Point{
    int x;
    int y;
};
struct TriPoint{
    int x;
    int y;
    int z;
};

double dist(position<double> a,position<double> b);
position<double> unit_Vector(position<double> a);
double CalcAngle(position<double> a,position<double> b, position<double> c);
double CalcCos(position<double> a,position<double> b, position<double> c);
TriPoint MakeTriangle(Point a,Point b);

#define CAPSULE 6
#define NANO 1

int main(int argc, char **argv)
{
	if(argc<3)
	{
		cout << "Usage: " << argv[0] << " radius AR\n";
		return 0;
	}

	
	
	
	//NP parameters
	int para_radius=atoi(argv[1]);
	double para_AspRat=atof(argv[2]);
	

    std::vector<position<double>> cylinderVertices;
    std::vector<position<double>> RcylinderVertices;
    std::vector<position<double>> LcylinderVertices;
    std::vector<position<double>> temp_LcylinderVertices;
    std::vector<position<double>> capsuleVertices;
    std::vector<position<double>> midPoint;
    std::vector<double> EdgeLength;
    std::vector<Point> EdgeList;
    std::vector<TriPoint> TriangleList;
    
    
    //Parameters
    double radius=10;					//Capsule is always constructed first with this radius
    double rad_scale=para_radius/radius;			//Actual radius=10*rad_scale
    double bondlength=1.5;
    double max_bondlength=1.8;				
    int AspRat=ceil(4/3*radius*(para_AspRat-1));	//No. of points inserted between the hemispheres
    double length=AspRat*bondlength;
    int Ntes=3;
    position<double> sphereCenter;
    sphereCenter.x = 0;
    sphereCenter.y = 0;
    sphereCenter.z = 0; 

    tesselaSphere<double> FirstSphere(sphereCenter, radius, Ntes);	//Generates the Sphere
    
    auto sphereVertex = FirstSphere.PrintSphereVertex();
        
    //SphereEnds for the Capsule
    for(int i = 0; i < sphereVertex.size(); i++)
    {
        position<double> temp;
	if (sphereVertex[i].x<sphereCenter.x){
            capsuleVertices.push_back(sphereVertex[i]);
        }
	if (sphereVertex[i].x>sphereCenter.x){
           temp.x = sphereVertex[i].x+length;
           temp.y = sphereVertex[i].y;
           temp.z = sphereVertex[i].z;
	   capsuleVertices.push_back(temp);
        }
     }
    
     
    //Left Ring for the cylinder
    std::vector<position<double>> Temp;
    for(int i = 0; i < sphereVertex.size(); i++)
    {
	if (sphereVertex[i].x==sphereCenter.x){
	    Temp.push_back(sphereVertex[i]);	    
    	}    
    }
    double index[Temp.size()][2];
    for(int i = 0; i < Temp.size(); i++)
    {
	double angle=atan2(Temp[i].y,Temp[i].z);
	index[i][0]=i;
	index[i][1]=angle;       
    }
    
    
    //Sort the points in clockwise direction
    double temp_f;
    double temp_s;
    for(int i = 0; i < Temp.size(); i++)
    {
    	for(int j = (i+1); j < Temp.size(); j++)
    	    {if (index[i][1]>index[j][1]){
	    	temp_f=index[j][0];
		temp_s=index[j][1];
		index[j][0]=index[i][0];
		index[j][1]=index[i][1];
		index[i][0]=temp_f;
		index[i][1]=temp_s;
		}
	     }
    }
    for(int i = 0; i < Temp.size(); i++)
    {	
	temp_LcylinderVertices.push_back(Temp[int(index[i][0])]);
    }
    
    //calculate the arc lengths for the points on the left ring
    double arc_length[temp_LcylinderVertices.size()];
    double min=100000;
    double max=0;
    for(int l = 0; l <temp_LcylinderVertices.size(); l++)
    {    if (l==temp_LcylinderVertices.size()-1){arc_length[l]=radius*abs(index[l][1]-((index[0][1])+2*M_PI));}
    	 else {arc_length[l]=radius*abs(index[l][1]-index[l+1][1]);}
	 if (arc_length[l]<min){min=arc_length[l];}
	 if (arc_length[l]>max){max=arc_length[l];}
    }
    
    //adjust the small arclengths on the left ring
    temp_f=0; //To calculate the arclength of actual points with small distance 
    int counter=0;  //Number of points with small distances to be adjusted
    int avg_counter=0;
    double a=0;		//Preferred distance between two points
    
    for(int l = 0; l <temp_LcylinderVertices.size(); l++){
    	if (arc_length[l]<=(min+0.1*min)){
		temp_f+=arc_length[l];
		counter=counter+1;
	}
	else { a+=arc_length[l];
		avg_counter+=1;}
    }
    
    
    a=max;
    int new_points=ceil(temp_f/(2*a)); //for each time
    double newLen=new_points*a;		//length of each time
    double epsilon=newLen-(temp_f/2.0);	//extra length each time
    
    
    //Final Left Ring for the cylinder
    double theta_prime=(a-epsilon/double(new_points))/radius;
    
    for(int l = 0; l <temp_LcylinderVertices.size(); l++){
    	if(l==(temp_LcylinderVertices.size()-1)){
	    if (arc_length[l]>(min+0.1*min))
	    {   LcylinderVertices.push_back(temp_LcylinderVertices[0]);
	    }
	    if (arc_length[l]<=(min+0.1*min))
	    {   
	        position<double> prev_temp;
		position<double> temp;
		prev_temp.x=0;
		prev_temp.y=temp_LcylinderVertices[l].y;
	        prev_temp.z=temp_LcylinderVertices[l].z;
		
		for (int i=0; i<new_points;i++){
		    temp.x=0;
		    temp.y=prev_temp.y*cos(-theta_prime)-prev_temp.z*sin(-theta_prime);
	            temp.z=prev_temp.y*sin(-theta_prime)+prev_temp.z*cos(-theta_prime);
		    LcylinderVertices.push_back(temp);
		    prev_temp=temp;
	        }
	        l=l+int(counter/2)-1;
	    }
	}
	else{
	    if (arc_length[l]>(min+0.1*min))
	    {   LcylinderVertices.push_back(temp_LcylinderVertices[l+1]);
	        temp_f=0;
	    }
	    if (arc_length[l]<=(min+0.1*min))
	    {   
	        position<double> prev_temp;
		position<double> temp;
		prev_temp.x=0;
		prev_temp.y=temp_LcylinderVertices[l].y;
	        prev_temp.z=temp_LcylinderVertices[l].z;
		
		for (int i=0; i<new_points;i++){
		    temp.x=0;
		    temp.y=prev_temp.y*cos(-theta_prime)-prev_temp.z*sin(-theta_prime);
	            temp.z=prev_temp.y*sin(-theta_prime)+prev_temp.z*cos(-theta_prime);
		    LcylinderVertices.push_back(temp);
		    prev_temp=temp;
	        }
	        l=l+ceil(counter/2)-1;
	    }
	}	
    }
    
    
     //Image of Right Ring on the cylinder
     for(int i = 0; i < LcylinderVertices.size(); i++)
     { 	 
         position<double> temp;
	 temp.x = LcylinderVertices[i].x+length;
         temp.y = LcylinderVertices[i].y;
         temp.z = LcylinderVertices[i].z;
	 RcylinderVertices.push_back(temp);
     }
     
                
    //Joining the points on the Ring
    for(int i = 0; i < LcylinderVertices.size(); i++)
    {
         position<double> nextP;
	 position<double> currentP;
	 currentP=LcylinderVertices[i];
	
	 while (currentP.x<=RcylinderVertices[i].x)
	 {
	 	
		cylinderVertices.push_back(currentP);
		capsuleVertices.push_back(currentP);//for the final vertices of capsule
	        nextP.x=currentP.x+bondlength;
		nextP.y=LcylinderVertices[i].y;
	 	nextP.z=LcylinderVertices[i].z;
		currentP=nextP;
		
         }
	 
     }
     
     //Midpoint on the cylinder between the lines
     
     for(int i = 0; i < cylinderVertices.size(); i++)
     {
	 position<double> mid;
	 position<double> vect;
	 int j= i%(AspRat+1); //there is AspRat+1 points in one straight line
	 int k=int(i/(AspRat+1));
	 if (j<AspRat)
	 {
		
		if (k<(LcylinderVertices.size()-1)){
		
		    mid.x=(cylinderVertices[i].x+cylinderVertices[i+AspRat+2].x)/2.0;
		    mid.y=(cylinderVertices[i].y+cylinderVertices[i+AspRat+2].y)/2.0;
	            mid.z=(cylinderVertices[i].z+cylinderVertices[i+AspRat+2].z)/2.0;
		    vect.x=0;
		    vect.y=mid.y;
		    vect.z=mid.z;
		    vect=unit_Vector(vect);
		    mid.y=radius*vect.y;
		    mid.z=radius*vect.z;
		    midPoint.push_back(mid);
		    capsuleVertices.push_back(mid);//for the final vertices of capsule
       		}
		if (k==(LcylinderVertices.size()-1)){
		    mid.x=(cylinderVertices[i].x+cylinderVertices[j+1].x)/2.0;
		    mid.y=(cylinderVertices[i].y+cylinderVertices[j+1].y)/2.0;
	            mid.z=(cylinderVertices[i].z+cylinderVertices[j+1].z)/2.0;
	 	    vect.x=0;
		    vect.y=mid.y;
		    vect.z=mid.z;
		    vect=unit_Vector(vect);
		    mid.y=radius*vect.y;
		    mid.z=radius*vect.z;
		    midPoint.push_back(mid);
		    capsuleVertices.push_back(mid);//for the final vertices of capsule
		}
         }
    }
    
    //join the cylinder and the sphere
    int k_1,k_2;
    for(int i = 0; i < LcylinderVertices.size(); i++)
     {
	counter=0;
	int k_1index[2];
	int k_2index[2];
	k_1=0;
	k_2=0;
	//joining point between cylinder and sphere to the left
	for  (int l = 0; l < sphereVertex.size(); l++)
	{   if (sphereVertex[l].x<sphereCenter.x){
		double distance=dist(sphereVertex[l],LcylinderVertices[i]);
		if (distance<1.8)
		{	counter+=1;
			if (k_1<2){
				k_1index[k_1]=l;
				k_1+=1;
			}
			
		}
		distance=dist(sphereVertex[l],LcylinderVertices[i+1]);
		if (distance<1.8)
		{	counter+=1;
			if (k_2<2){
				k_2index[k_2]=l;
				k_2+=1;
			}
		}
	      }
	      if (counter==5 || counter==6){break;}
	}
	min=100000; 		//just a arbitrary number
	if (counter ==4)
	{	position<double> temp;
		position<double> vect;
		//calculate two points on the sphere to join
		double dist1=dist(sphereVertex[k_1index[0]],sphereVertex[k_2index[0]]);
		double dist2=dist(sphereVertex[k_1index[0]],sphereVertex[k_2index[1]]);
		double dist3=dist(sphereVertex[k_1index[1]],sphereVertex[k_2index[0]]);
		double dist4=dist(sphereVertex[k_1index[1]],sphereVertex[k_2index[1]]);
		if (dist1<dist2 && dist1<dist3 && dist1<dist4){k_1=k_1index[0];k_2=k_2index[0];}
		if (dist2<dist1 && dist2<dist3 && dist2<dist4){k_1=k_1index[0];k_2=k_2index[1];}
		if (dist3<dist2 && dist3<dist1 && dist3<dist4){k_1=k_1index[1];k_2=k_2index[0];}
		if (dist4<dist1 && dist4<dist2 && dist4<dist3){k_1=k_1index[1];k_2=k_2index[1];}
		int j= i*(AspRat+1);
		temp.x=(sphereVertex[k_1].x+sphereVertex[k_2].x+LcylinderVertices[i].x+LcylinderVertices[i+1].x+midPoint[j-i].x)/5;
		temp.y=(sphereVertex[k_1].y+sphereVertex[k_2].y+LcylinderVertices[i].y+LcylinderVertices[i+1].y+midPoint[j-i].y)/5;
		temp.z=(sphereVertex[k_1].z+sphereVertex[k_2].z+LcylinderVertices[i].z+LcylinderVertices[i+1].z+midPoint[j-i].z)/5;
		
		vect.x=temp.x-sphereCenter.x;
		vect.y=temp.y-sphereCenter.y;
		vect.z=temp.z-sphereCenter.z;
		vect=unit_Vector(vect);
		temp.x=radius*vect.x;
		temp.y=radius*vect.y;
		temp.z=radius*vect.z;
		capsuleVertices.push_back(temp); //join the points on the left
		temp.x=sphereCenter.x+length+abs(sphereCenter.x-temp.x);
		capsuleVertices.push_back(temp);  //join the points on the right
			
      }
   }
   
   //EdgeMap bool to store the Edge information
   int n=capsuleVertices.size();
   bool E_map[n][n];
   for (int i = 0; i < n; i++){
   	for (int j = 0; j < n; j++){
		E_map[i][j] = false;
		if (i!=j)
		{	double distance=dist(capsuleVertices[i],capsuleVertices[j]);
			if (distance<max_bondlength){E_map[i][j]=true;}
		}
	}
   }
   
   //Temporary store the indices of the Edges, and the Edgelengths
   Point EdgePoint;
   std::vector<Point> Temp_List;
   for (int i = 0; i <n; i++){
   	for (int j = 0; j < i; j++){
   		if (E_map[i][j])
		{	
			
			EdgePoint.x=i;
			EdgePoint.y=j;
			Temp_List.push_back(EdgePoint);
			EdgeLength.push_back(dist(capsuleVertices[i],capsuleVertices[j]));
		}	
	}
   }
   		
   
   
   //Store the Edges after rearranging it and then store the Triangles
   //Temp triangle to store all the edges connected to each point, and then store them as triangles
   std::vector<Point> TempTriangle;
   for (int j=0;j<n;j++){
   	for (int i=0; i<Temp_List.size();i++){	
   		if (j==Temp_List[i].x){
			EdgeList.push_back(Temp_List[i]);
			TempTriangle.push_back(Temp_List[i]);
		}
		if (j==Temp_List[i].y){
			EdgePoint.x=Temp_List[i].y;
			EdgePoint.y=Temp_List[i].x;
			EdgeList.push_back(EdgePoint);
			TempTriangle.push_back(EdgePoint);
		}
   	}
	
	int index;
	double min=10000;
	int l=0;
	EdgePoint=TempTriangle[0];
	
	
	while (!TempTriangle.empty()){
		for (int k=1; k<TempTriangle.size();k++){
			double theta=CalcAngle(capsuleVertices[TempTriangle[0].x],capsuleVertices[TempTriangle[0].y], capsuleVertices[TempTriangle[k].y]);
		    	if (theta<min){min=theta;index=k;}
		}
		TriangleList.push_back(MakeTriangle(TempTriangle[0],TempTriangle[index]));
		TempTriangle.erase (TempTriangle.begin());
		if (TempTriangle.size()==1)
		{	TriangleList.push_back(MakeTriangle(TempTriangle[0],EdgePoint));
			TempTriangle.clear();
		}
	}
   }
   
   	//change the radius of the particle
	for(int i=0;i<capsuleVertices.size();i++)
	{
		capsuleVertices[i].x=capsuleVertices[i].x*rad_scale;
		capsuleVertices[i].y=capsuleVertices[i].y*rad_scale;
		capsuleVertices[i].z=capsuleVertices[i].z*rad_scale;
	}
	
	
	ofstream fout;
    ostringstream outstr;
    outstr <<"Capsule_rad_"<<para_radius<<"_AR_"<<para_AspRat<<"_vertices.xyz";
    string ofilename=outstr.str();
    fout.open(ofilename.c_str(), std::ios::out);
    if(fout.is_open()){
         fout <<capsuleVertices.size()<<endl;
         fout <<"test" <<endl;
         for(int i = 0; i < capsuleVertices.size(); i++)
   	 {
		fout<<0<<" "<<capsuleVertices[i].x<<" "<<capsuleVertices[i].y<<"  "<<capsuleVertices[i].z<<" "<< std::endl;
   	 }
    }
    fout.close();

	return 0;
}

double dist(position<double> a,position<double> b)
{       double dx= a.x-b.x;
	double dy= a.y-b.y;
	double dz= a.z-b.z;
	double distance=dx*dx+dy*dy+dz*dz;
	return powf(distance,0.5);
}

position<double> unit_Vector(position<double> a)
{
	double m=a.x*a.x+a.y*a.y+a.z*a.z;
	m=powf(m,0.5);
	a.x/=m;
	a.y/=m;
	a.z/=m;
	return a;
}

double CalcAngle(position<double> a,position<double> b, position<double> c)
{
	double dx1= b.x-a.x;
	double dy1= b.y-a.y;
	double dz1= b.z-a.z;
	double distance1=dx1*dx1+dy1*dy1+dz1*dz1;
	distance1=powf(distance1,0.5);
		
	double dx2= c.x-a.x;
	double dy2= c.y-a.y;
	double dz2= c.z-a.z;
	double distance2=dx2*dx2+dy2*dy2+dz2*dz2;
	distance2=powf(distance2,0.5);
	
	double cos_theta=(dx1*dx2+dy1*dy2+dz1*dz2)/(distance1*distance2);
	return acos(cos_theta);
	
}

double CalcCos(position<double> a,position<double> b, position<double> c)
{
	double dx1= b.x-a.x;
	double dy1= b.y-a.y;
	double dz1= b.z-a.z;
	double distance1=dx1*dx1+dy1*dy1+dz1*dz1;
	distance1=powf(distance1,0.5);
		
	double dx2= c.x-a.x;
	double dy2= c.y-a.y;
	double dz2= c.z-a.z;
	double distance2=dx2*dx2+dy2*dy2+dz2*dz2;
	distance2=powf(distance2,0.5);
	
	double cos_theta=(dx1*dx2+dy1*dy2+dz1*dz2)/(distance1*distance2);
	return cos_theta;
	
}

TriPoint MakeTriangle(Point a,Point b)
{
	TriPoint temp;
	temp.x=a.x;
	temp.y=a.y;
	temp.z=b.y;
	return temp;
}
