#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <stdint.h>

void cross(double a[3],double b[3],double *c);
double dot(double a[3],double b[3]);

int main(int argc, char *argv[])
{
	FILE *inputfpt, *outputfpt;	
	
	unsigned int vertex = 0;
	unsigned int face = 0;
	unsigned int vert;
	unsigned int fac;
	unsigned char *end;
	unsigned char *fptr;
	unsigned char temp;	
	unsigned char extension[3];
	
	double camera[] = {1.0, 0.0, 0.0};
	double up[] = { 0.0, 0.0, 1.0};
	double xt, yt, zt;
	double xr, yr, zr;
	
	float E = 0;	
	float **verticies;
	float current[3];	
	float min[] = {0.0, 0.0, 0.0};
	float max[] = {0.0, 0.0, 0.0};
	float center[] = {0.0, 0.0, 0.0};	
	
	int **faces;
	int start = 0;
	int size = 0;
	int i,j,k;
	int tracker = 0;
	int tx, ty, tz;
	int curr[3];	
	

	//===================================================================//
	//                     open input and output files                   //	
	//===================================================================//
	if(argc != 5)
	{
		printf("Error: Invalid number of arguments\n");
		exit(0);
	}
		
	//get degrees from command line
	tx = atoi(argv[2]);
	ty = atoi(argv[3]);
	tz = atoi(argv[4]);
	xt = (double)tx;
	yt = (double)ty;
	zt = (double)tz;

	//convert degrees to radians
	xr = xt*M_PI/180;
	yr = yt*M_PI/180;
	zr = zt*M_PI/180;
	printf("x: %fl\n", xt);
	printf("y: %fl\n", yt);
	printf("z: %fl\n", zt);

	fptr = argv[1];
	inputfpt = fopen(fptr, "r");
	if(inputfpt == NULL)
	{
		printf("Unable to open %s\n",fptr);
		exit(0);
	}	
	end = fptr + strlen(fptr);
	while(end > fptr && *end != '.')
	{
		--end;
	}
	if(end > fptr)
	{
		*end = '\0';
	}
	strcpy(extension,".ppm");
	strcat(fptr,extension);
	printf("Output file: %s\n",fptr);
	printf("\n");
	outputfpt = fopen(fptr, "wb");
	if(outputfpt == NULL)
	{
		printf("Unable to open %s\n",fptr);
		exit(0);
	}

	//===================================================================//
	//                        get # of verticies                         //	
	//===================================================================//
	while((temp = fgetc(inputfpt)) != 'x');
	temp = fgetc(inputfpt);
	while((temp = fgetc(inputfpt)) != '\n')
	{
		size++;
	}
	size -= 1;
	rewind(inputfpt);
	while((temp = fgetc(inputfpt)) != 'x');
	temp = fgetc(inputfpt);
	for(i = 0; i < size; i++)
	{
		temp = fgetc(inputfpt);
		vert = temp - '0';
		vertex = vertex + (vert * pow(10, (size-i-1)));	
	}
	printf("vertex: %d\n",vertex);
	
	//===================================================================//
	//                          get # of faces                           //	
	//===================================================================//
	rewind(inputfpt);
	while(tracker < 6)
	{
		if((temp = fgetc(inputfpt)) == '\n')
		{
			tracker++;
		}
	}
	for(i = 0; i < 13; i++)
	{
		temp = fgetc(inputfpt);
	}
	size = 0;
	while((temp = fgetc(inputfpt)) != '\n')
	{
		size++;
	}
	size -= 1;
	rewind(inputfpt);
	tracker = 0;
	while(tracker < 6)
	{
		if((temp = fgetc(inputfpt)) == '\n')
		{
			tracker++;
		}
	}
	for(i = 0; i < 13; i++)
	{
		temp = fgetc(inputfpt);
	}
	for(i = 0; i < size; i++)
	{
		temp = fgetc(inputfpt);
		fac = temp - '0';
		face = face + (fac * pow(10, (size-i-1)));		
	}
	printf("face: %d\n",face);

	//===================================================================//
	//                     read verticies and faces                      //	
	//===================================================================//
	verticies = (float **)calloc(vertex, sizeof(float *));
	for(i = 0; i < vertex; i++)
	{
		verticies[i] = (float *)calloc(3, sizeof(float));
	}
	faces = (int **)calloc(face, sizeof(int *));
	for(i = 0; i < face; i++)
	{
		faces[i] = (int *)calloc(3, sizeof(int));
	}
	
	rewind(inputfpt);
	tracker = 0;
	while(tracker < 9)
	{
		if((temp = fgetc(inputfpt)) == '\n')
		{
			tracker++;
		}
	}

	for(i = 0; i < vertex; i++)
	{
		fscanf(inputfpt, "%f", &current[0]);
		fscanf(inputfpt, "%f", &current[1]);
		fscanf(inputfpt, "%f", &current[2]);
		verticies[i][0] = current[0];
		verticies[i][1] = current[1];
		verticies[i][2] = current[2];
		
		if(start == 0)
		{
			max[0] = current[0];
			max[1] = current[1];
			max[2] = current[2];
			min[0] = current[0];
			min[1] = current[1];
			min[2] = current[2];
			start = 1;
		}		

		//find max
		if(current[0] > max[0])
		{
			max[0] = current[0];
		}
		if(current[1] > max[1])
		{
			max[1] = current[1];
		}
		if(current[2] > max[2])
		{
			max[2] = current[2];
		}

		//find min
		if(current[0] < min[0])
		{
			min[0] = current[0];
		}
		if(current[1] < min[1])
		{
			min[1] = current[1];
		}
		if(current[2] < min[2])
		{
			min[2] = current[2];
		}
	}

	//calculate center
	center[0] = (max[0] + min[0]) / 2;
	center[1] = (max[1] + min[1]) / 2;
	center[2] = (max[2] + min[2]) / 2;
	
	//find extent, E
	E = max[0] - min[0];
	if((max[1] - min[1]) > E)
	{
		E = (max[1] - min[1]);
	}
	if((max[2] - min[2]) > E)
	{
		E = (max[2] - min[2]);
	}

	for(i = 0; i < face; i++)
	{
		fscanf(inputfpt, "%d", &curr[0]);
		fscanf(inputfpt, "%d", &curr[0]);
		fscanf(inputfpt, "%d", &curr[1]);
		fscanf(inputfpt, "%d", &curr[2]);
		faces[i][0] = curr[0];
		faces[i][1] = curr[1];
		faces[i][2] = curr[2];
	}

	//===================================================================//
	//                      Camera and Up rotations                      //	
	//===================================================================//
	double left[3];
	double a;
	double right[3];
	double top[3];
	double bottom[3];
	double topleft[3];
	double diff[3];
	double buff[3];
	double buff2[3];
	double rx[3][3] = {{1,0,0},{0,cos(xr),sin(xr)},{0,-sin(xr),cos(xr)}};
	double ry[3][3] = {{cos(yr),0,-sin(yr)},{0,1,0},{sin(yr),0,cos(yr)}};
	double rz[3][3] = {{cos(zr),sin(zr),0},{-sin(zr),cos(zr),0},{0,0,1}};

	//rotations for camera and up	
	for(i = 0; i < 3; i++)
	{
		buff[i] = camera[0]*rx[0][i] + camera[1]*rx[1][i] + camera[2]*rx[2][i];
		buff2[i] = up[0]*rx[0][i] + up[1]*rx[1][i] + up[2]*rx[2][i];
	}
	for(i = 0; i < 3; i++)
	{
		camera[i] = buff[i];
		up[i] = buff2[i];
	}
	for(i = 0; i < 3; i++)
	{
		buff[i] = camera[0]*ry[0][i] + camera[1]*ry[1][i] + camera[2]*ry[2][i];
		buff2[i] = up[0]*ry[0][i] + up[1]*ry[1][i] + up[2]*ry[2][i];
	}
	for(i = 0; i < 3; i++)
	{
		camera[i] = buff[i];
		up[i] = buff2[i];
	}
	for(i = 0; i < 3; i++)
	{
		buff[i] = camera[0]*rz[0][i] + camera[1]*rz[1][i] + camera[2]*rz[2][i];
		buff2[i] = up[0]*rz[0][i] + up[1]*rz[1][i] + up[2]*rz[2][i];
	}
	for(i = 0; i < 3; i++)
	{
		camera[i] = buff[i];
		up[i] = buff2[i];
	}

	for(i = 0; i < 3; i++)
	{
		camera[i] = 1.5*E*camera[i] + center[i];
	}
	
	//===================================================================//
	//                             Bound Box                             //	
	//===================================================================//
	for(i = 0; i < 3; i++)
	{
		diff[i] = center[i]-camera[i];
	}
	left[0] = up[1]*diff[2] - up[2]*diff[1];	
	left[1] = up[2]*diff[0] - up[0]*diff[2];	
	left[2] = up[0]*diff[1] - up[1]*diff[0];	
	
	a = sqrt( pow(left[0],2) + pow(left[1],2) + pow(left[2],2) );

	for(i = 0; i < 3; i++)
	{
		left[i] = (E/(2*a))*left[i]+center[i];
	}	

	for(i = 0; i < 3; i++)
	{
		diff[i] =  center[i] - camera[i]; 
	}
	right[0] = diff[1]*up[2] - diff[2]*up[1];
	right[1] = diff[2]*up[0] - diff[0]*up[2];
	right[2] = diff[0]*up[1] - diff[1]*up[0];

	for(i = 0; i < 3; i++)
	{
		right[i] = (E/(2*a))*right[i] + center[i];
	}

	for(i = 0; i < 3; i++)
	{
		top[i] = (E/2)*up[i]+center[i];
	}

	for(i = 0; i < 3; i++)
	{
		bottom[i] = (-E/2)*up[i] + center[i];
	}	

	for(i = 0; i < 3; i++)
	{
		topleft[i] = (E/2)*up[i] + left[i];
	}

	printf("\ncamera: %.2f %.2f %.2f \n", camera[0], camera[1], camera[2]);
	printf("up: %.2f %.2f %.2f \n", up[0], up[1], up[2]);	
	printf("Extent: %.4f\n", E);
	printf("\n");
	printf("Left: %.2f %.2f %.2f \n", left[0], left[1], left[2]);
	printf("Right: %.2f %.2f %.2f \n", right[0], right[1], right[2]);
	printf("\n");	
	printf("Top: %.2f %.2f %.2f \n", top[0], top[1], top[2]);
	printf("Bottom: %.2f %.2f %.2f \n", bottom[0], bottom[1], bottom[2]);
	printf("\n");
	printf("TopLeft: %.2f %.2f %.2f \n", topleft[0], topleft[1], topleft[2]);
	
	//===================================================================//
	//                          Calculations                             //	
	//===================================================================//
	int COLS = 256, ROWS = 256;
	int c, r;
	int f0, f1, f2;

	double zbuffer = 999999;
	unsigned char color[ROWS][COLS];
	
	double *plane;
	double capD;
	double n, d;
	double cmul, rmul;

	double dot1, dot2, dot3;
	double *cross1, *cross2;
	double intersect[3];	
	double image[3];
	double *diff1, *diff2, *diff3, *diff4;			
	double v0[3], v1[3], v2[3];
	
	plane = (double *)calloc(3, sizeof(double));
	cross1 = (double *)calloc(3, sizeof(double));
	cross2 = (double *)calloc(3, sizeof(double));
	diff1 = (double *)calloc(3, sizeof(double));
	diff2 = (double *)calloc(3, sizeof(double));
	diff3 = (double *)calloc(3, sizeof(double));
	diff4 = (double *)calloc(3, sizeof(double));

	for(r = 0; r < ROWS; r++)
	{
		for(c = 0; c < COLS; c++)
		{
			color[r][c] = 0;
			zbuffer = 999999;
			cmul = (double)c/(COLS-1);
			rmul = (double)r/(ROWS-1);
			
			for(i = 0; i < 3; i++)
			{
				diff1[i] = right[i]-left[i];
				diff2[i] = bottom[i]-top[i];
				image[i] = topleft[i] + cmul*diff1[i] + rmul*diff2[i];
			}		

			//for each triangle with coordinates v0 v1 v2
			for(i = 0; i < face; i++)
			{
				f0 = faces[i][0];
				f1 = faces[i][1];
				f2 = faces[i][2];

				//for each triangle, get its verticies 
				for(j = 0; j < 3; j++)
				{
					v0[j] = verticies[f0][j];
					v1[j] = verticies[f1][j];
					v2[j] = verticies[f2][j];
				}

				for(j = 0; j < 3; j++)
				{
					diff1[j] = v1[j] - v0[j];
					diff2[j] = v2[j] - v0[j];
				}			
				cross(diff1,diff2,plane);

				for(j = 0; j < 3; j++)
				{
					buff[j] = -plane[j];
				}
				capD = dot(buff,v0);
				
				n = dot(buff,camera) - capD;

				for(j = 0; j < 3; j++)
				{
					diff[j] = image[j] - camera[j];
				}
				d = dot(plane,diff);

				if(d != 0)
				{
					for(j = 0; j < 3; j++)
					{
						intersect[j] = camera[j] + (n/d) * (image[j] - camera[j]);	
					}
					
					for(j = 0; j < 3; j++)
					{
						diff1[j] = v2[j] - v0[j];
						diff2[j] = v1[j] - v0[j];
						diff3[j] = intersect[j] - v0[j];	
					}

					cross1[0] = diff1[1]*diff2[2] - diff1[2]*diff2[1];
					cross1[1] = diff1[2]*diff2[0] - diff1[0]*diff2[2];
					cross1[2] = diff1[0]*diff2[1] - diff1[1]*diff2[0];
					cross2[0] = diff3[1]*diff2[2] - diff3[2]*diff2[1];
					cross2[1] = diff3[2]*diff2[0] - diff3[0]*diff2[2];
					cross2[2] = diff3[0]*diff2[1] - diff3[1]*diff2[0];
					dot1 = dot(cross1,cross2);				 

					for(j = 0; j < 3; j++)
					{
						diff1[j] = v0[j] - v1[j];
						diff2[j] = v2[j] - v1[j];
						diff3[j] = intersect[j] - v1[j];	
					}
					cross1[0] = diff1[1]*diff2[2] - diff1[2]*diff2[1];
					cross1[1] = diff1[2]*diff2[0] - diff1[0]*diff2[2];
					cross1[2] = diff1[0]*diff2[1] - diff1[1]*diff2[0];
					cross2[0] = diff3[1]*diff2[2] - diff3[2]*diff2[1];
					cross2[1] = diff3[2]*diff2[0] - diff3[0]*diff2[2];
					cross2[2] = diff3[0]*diff2[1] - diff3[1]*diff2[0];
					dot2 = dot(cross1,cross2);				 

					for(j = 0; j < 3; j++)
					{
						diff1[j] = v1[j] - v2[j];
						diff2[j] = v0[j] - v2[j];
						diff3[j] = intersect[j] - v2[j];	
					}

					cross1[0] = diff1[1]*diff2[2] - diff1[2]*diff2[1];
					cross1[1] = diff1[2]*diff2[0] - diff1[0]*diff2[2];
					cross1[2] = diff1[0]*diff2[1] - diff1[1]*diff2[0];
					cross2[0] = diff3[1]*diff2[2] - diff3[2]*diff2[1];
					cross2[1] = diff3[2]*diff2[0] - diff3[0]*diff2[2];
					cross2[2] = diff3[0]*diff2[1] - diff3[1]*diff2[0];
					dot3 = dot(cross1,cross2);	
					
					if(((dot1 >= 0) && (dot2 >= 0) && (dot3 >= 0) ) && ((n/d) <= zbuffer) )
					{
						zbuffer = n/d;
						color[r][c] = 155 + (i % 100);			
					}
				
			 
				}
	
			}
			

		}
	}
	//===================================================================//
	//                         Create PPM image                          //	
	//===================================================================//
	
	(void) fprintf(outputfpt, "P5\n%d %d\n255\n",256,256);
	(void) fwrite(&color,sizeof(unsigned char), ROWS*COLS,outputfpt);	

	//free and close
	for(i = 0; i < vertex; i++)
	{
		free(verticies[i]);
	}
	free(verticies);
	for(i = 0; i < face; i++)
	{
		free(faces[i]);
	}
	free(faces);
	free(diff1);
	free(diff2);
	free(diff3);
	free(diff4);
	free(cross1);
	free(cross2);
	free(plane);
	fclose(inputfpt);
	fclose(outputfpt);
	return 0;
}

void cross(double a[3],double b[3],double *c)
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];

}

double dot(double a[3],double b[3])
{
	int i;
	double c = 0;
	for(i = 0; i < 3; i++)
	{
		c += a[i]*b[i];	
	}

	return c;
}


