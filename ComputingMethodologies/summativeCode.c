// (C) 2017 Tobias Weinzierl edits by Joseph Cass

#include <fstream>
#include <sstream>
#include <iostream>
#include <omp.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>


double t = 0;
double tFinal = 0;

int NumberOfBodies = 0;

double** x;
double** v;
double*  mass;


void setUp(int argc, char** argv) {


	if (argc == 3) { //random bodies input
		tFinal = std::stof(argv[1]);
		NumberOfBodies = std::stof(argv[2]);
		x = new double*[NumberOfBodies]; //pointer
		v = new double*[NumberOfBodies]; //pointer
		mass = new double[NumberOfBodies]; //pointer

		int minVel = -10;
		int maxVel = 10;
		int minPos = -50;
		int maxPos = 50;
		int minMass = 1;
		int maxMass = 1;
		for (int i = 0; i < NumberOfBodies; i++)
		{
			x[i] = new double[3];
			v[i] = new double[3];

			x[i][0] = rand() % (maxPos + 1 - minPos) + minPos;
			x[i][1] = rand() % (maxPos + 1 - minPos) + minPos;
			x[i][2] = rand() % (maxPos + 1 - minPos) + minPos;
			v[i][0] = rand() % (maxVel + 1 - minVel) + minVel;
			v[i][1] = rand() % (maxVel + 1 - minVel) + minVel;
			v[i][2] = rand() % (maxVel + 1 - minVel) + minVel;
			mass[i] = rand() % (maxMass + 1 - minMass) + minMass;
		}

		std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

	}

	else {
		NumberOfBodies = (argc - 2) / 7;
		x = new double*[NumberOfBodies]; //pointer
		v = new double*[NumberOfBodies]; //pointer
		mass = new double[NumberOfBodies]; //pointer
		int readArgument = 1;


		tFinal = std::stof(argv[readArgument]); readArgument++;
		printf("tFinal: %lf \n", tFinal);

		for (int i = 0; i < NumberOfBodies; i++) {

			x[i] = new double[3];
			v[i] = new double[3];

			x[i][0] = std::stof(argv[readArgument]); readArgument++;
			x[i][1] = std::stof(argv[readArgument]); readArgument++;
			x[i][2] = std::stof(argv[readArgument]); readArgument++;

			v[i][0] = std::stof(argv[readArgument]); readArgument++;
			v[i][1] = std::stof(argv[readArgument]); readArgument++;
			v[i][2] = std::stof(argv[readArgument]); readArgument++;

			mass[i] = std::stof(argv[readArgument]); readArgument++;

			if (mass[i] <= 0.0) {
				std::cerr << "invalid mass for body " << i << std::endl;
				exit(-2);
			}
		}

		std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
	}
}


std::ofstream videoFile;


void openParaviewVideoFile() {
	videoFile.open("result.pvd");
	videoFile << "<?xml version=\"1.0\"?>" << std::endl
		<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
		<< "<Collection>";
}



void closeParaviewVideoFile() {
	videoFile << "</Collection>"
		<< "</VTKFile>" << std::endl;
}


/**
* The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
*/
void printParaviewSnapshot(int counter) {
	std::stringstream filename;
	filename << "result-" << counter << ".vtp";
	std::ofstream out(filename.str().c_str());
	out << "<VTKFile type=\"PolyData\" >" << std::endl
		<< "<PolyData>" << std::endl
		<< " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
		<< "  <Points>" << std::endl
		<< "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

	for (int i = 0; i < NumberOfBodies; i++) {
		out << x[i][0]
			<< " "
			<< x[i][1]
			<< " "
			<< x[i][2]
			<< " ";
	}

	out << "   </DataArray>" << std::endl
		<< "  </Points>" << std::endl
		<< " </Piece>" << std::endl
		<< "</PolyData>" << std::endl
		<< "</VTKFile>" << std::endl;

	videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}


void updateBody() {
	int NumberOfBodiesDeleted = 0;
	//*********************************************SETEP*********************************************
	double timeStepSize = 0.0001; //original
	//double timeStepSize = 1e-10;

	double forces[NumberOfBodies][3];				//create forces array

#pragma omp parallel
{
	#pragma omp for
	for (int body = 0; body < NumberOfBodies; body++)
	{
		forces[body][0] = 0;
		forces[body][1] = 0;
		forces[body][2] = 0;
	}

	int collisionOccured = 0;

	//***********************************************************************************************

	//*********************************************CALCULATE FORCES**********************************
	#pragma omp for
	for (int eachBody = 0; eachBody < NumberOfBodies; eachBody++) {						//for each body


		for (int eachOtherBody = 0; eachOtherBody < NumberOfBodies; eachOtherBody++) {	//compare against each other body
			if (eachBody != eachOtherBody) {											//Don't compare bodies against themselves
				double distance = sqrt(											//calculate distances from unchanged values
					(x[eachBody][0] - x[eachOtherBody][0]) * (x[eachBody][0] - x[eachOtherBody][0]) +
					(x[eachBody][1] - x[eachOtherBody][1]) * (x[eachBody][1] - x[eachOtherBody][1]) +
					(x[eachBody][2] - x[eachOtherBody][2]) * (x[eachBody][2] - x[eachOtherBody][2])
				);

				if (distance < 0.000000001  && eachBody != eachOtherBody) //0.000000001 for Windows, 0.000000001 + 0.00000000000001
				{
					printf("COLLISON \n");
					printf("Body %d:\n	Position:(%.30f,%lf,%lf)\n	Velocity:(%lf,%lf,%lf)\n	Mass:%lf\n", eachBody, x[eachBody][0], x[eachBody][1], x[eachBody][2], v[eachBody][0], v[eachBody][1], v[eachBody][2], mass[eachBody]);

					printf("Body %d:\n	Position:(%.30f,%lf,%lf)\n	Velocity:(%lf,%lf,%lf)\n	Mass:%lf\n", eachOtherBody, x[eachOtherBody][0], x[eachOtherBody][1], x[eachOtherBody][2], v[eachOtherBody][0], v[eachOtherBody][1], v[eachOtherBody][2], mass[eachOtherBody]);

					printf("\n");
					printf("Position of collision: (%.30f,%lf,%lf) \n", x[eachBody][0], x[eachBody][1], x[eachBody][2]);

					//merge boddies
					printf("\n");
					printf("Merged body: \n");
					//0: find momentums and average velocities:
					//m1*v1 + m2*v2 = (m1 + m2)*Vf
					double newVelocityX = (mass[eachBody] * v[eachBody][0] + mass[eachOtherBody] * v[eachOtherBody][0]) / (mass[eachBody] + mass[eachOtherBody]);
					double newVelocityY = (mass[eachBody] * v[eachBody][1] + mass[eachOtherBody] * v[eachOtherBody][1]) / (mass[eachBody] + mass[eachOtherBody]);
					double newVelocityZ = (mass[eachBody] * v[eachBody][2] + mass[eachOtherBody] * v[eachOtherBody][2]) / (mass[eachBody] + mass[eachOtherBody]);

					//2: adapt merged body
					v[eachBody][0] = newVelocityX;
					v[eachBody][1] = newVelocityY;
					v[eachBody][2] = newVelocityZ;

					//reset forces to zero?
					//try re-evalute force on body by deducting forces from subject bodies
					//forces[i][0] -= ((previousPosArray[j][0] - previousPosArray[i][0]) * mass[j] * mass[i] / distances[i][j] / distances[i][j] / distances[i][j]);
					forces[eachBody][0] = 0;
					forces[eachBody][1] = 0;
					forces[eachBody][2] = 0;

					mass[eachBody] += mass[eachOtherBody]; //increase mass
					printf("	Velocity: (%lf,%lf,%lf) \n	Mass:%lf\n", v[eachBody][0], v[eachBody][1], v[eachBody][2], mass[eachBody]);
					printf("\n \n");

					//3: delete one of the bodies
					// delete index j


					//{
						mass[eachOtherBody] = mass[NumberOfBodies - 1]; // copy next element left

						for (int direction = 0; direction < 3; direction++)
						{
							x[eachOtherBody][direction] = x[NumberOfBodies - 1][direction]; // copy next element left
							v[eachOtherBody][direction] = v[NumberOfBodies - 1][direction]; // copy next element left
							forces[eachOtherBody][direction] = forces[NumberOfBodies - 1][direction]; // copy next element left
						}
					//}

					//4: re-arrange compensate for change in array
					#pragma omp critical(critsect)
					{
					NumberOfBodies += -1;
					}
					
					distance = sqrt(											//calculate distances from unchanged values
						(x[eachBody][0] - x[eachOtherBody][0]) * (x[eachBody][0] - x[eachOtherBody][0]) +
						(x[eachBody][1] - x[eachOtherBody][1]) * (x[eachBody][1] - x[eachOtherBody][1]) +
						(x[eachBody][2] - x[eachOtherBody][2]) * (x[eachBody][2] - x[eachOtherBody][2])
					);

					printf("\nNumber of Bodies: %d \n", NumberOfBodies);

				}
				//distances[eachBody][eachOtherBody] = distance;

				forces[eachBody][0] += (x[eachOtherBody][0] - x[eachBody][0]) * mass[eachOtherBody] * mass[eachBody] / distance / distance / distance;
				forces[eachBody][1] += (x[eachOtherBody][1] - x[eachBody][1]) * mass[eachOtherBody] * mass[eachBody] / distance / distance / distance;
				forces[eachBody][2] += (x[eachOtherBody][2] - x[eachBody][2]) * mass[eachOtherBody] * mass[eachBody] / distance / distance / distance;
			}
		 } //end of j loop

		  //----EVALUATE TIME STEP----
		  //double eps = 0.00005; //original
		double eps = 0.005;
		//if input velocity is zero causes problems!

		while ((abs(timeStepSize * forces[eachBody][0] / mass[eachBody]) / v[eachBody][0] > eps && (v[eachBody][0] > 0.0000001 || v[eachBody][0] < -0.0000001))
			|| (abs(timeStepSize * forces[eachBody][1] / mass[eachBody]) / v[eachBody][1] > eps && (v[eachBody][1] > 0.0000001 || v[eachBody][1] < -0.0000001))
			|| (abs(timeStepSize * forces[eachBody][2] / mass[eachBody]) / v[eachBody][2] > eps && (v[eachBody][2] > 0.0000001 || v[eachBody][2] < -0.0000001)))
			{
			//printf("Decreasing timeStepSize \n");
			timeStepSize /= 2.0;
		}


		/*  Alternat time evaluation Doesn't consider change in acceleration
		if (sqrt(2 *eps / abs(forces[eachBody][0] / mass[eachBody])) < timeStepSize)
		{
		timeStepSize = sqrt(2 * eps / abs(forces[eachBody][0] / mass[eachBody]));
		//printf("reduced timestep");
		}
		if (sqrt(2 *eps / abs(forces[eachBody][1] / mass[eachBody])) < timeStepSize)
		{
		timeStepSize = sqrt(2 * eps / abs(forces[eachBody][1] / mass[eachBody]));
		//printf("reduced timestep");
		}
		if (sqrt(2 *eps / abs(forces[eachBody][2] / mass[eachBody])) < timeStepSize)
		{
		timeStepSize = sqrt(2 * eps / abs(forces[eachBody][2] / mass[eachBody]));
		//printf("reduced timestep");
		}
		*/
		//--------------------------


}//end of i loop


	//********************************UPDATE POSITIONS AND VELOCITY**********************************

	#pragma omp for
	for (int eachBody = 0; eachBody < NumberOfBodies; eachBody++)
	{
		x[eachBody][0] = x[eachBody][0] + timeStepSize * v[eachBody][0];
		x[eachBody][1] = x[eachBody][1] + timeStepSize * v[eachBody][1];
		x[eachBody][2] = x[eachBody][2] + timeStepSize * v[eachBody][2];

		v[eachBody][0] = v[eachBody][0] + timeStepSize * forces[eachBody][0] / mass[eachBody];
		v[eachBody][1] = v[eachBody][1] + timeStepSize * forces[eachBody][1] / mass[eachBody];
		v[eachBody][2] = v[eachBody][2] + timeStepSize * forces[eachBody][2] / mass[eachBody];
	}
	t += timeStepSize;
	//***********************************************************************************************
}
}



int main(int argc, char** argv) {
	omp_set_num_threads(1);
	if (argc == 1) {
		std::cerr << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
			<< std::endl
			<< "Examples:" << std::endl
			<< "100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
			<< "100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
			<< "100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
			<< std::endl
			<< "In this naive code, only the first body moves" << std::endl;

		return -1;
	}

	else if ((argc - 2) % 7 != 0 && argc != 3) {
		std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
		return -2;
	}

	setUp(argc, argv);

	openParaviewVideoFile();

	printParaviewSnapshot(0);

	int timeStepsSinceLastPlot = 0;
	const int plotEveryKthStep = 100;

	clock_t start, end;
	double cpu_time_used;

	//More efficient method? iterate for x calulations then snapshot rather than checking each time
	start = clock();
	while (t <= tFinal) {
		updateBody();
		timeStepsSinceLastPlot++;

		if (timeStepsSinceLastPlot%plotEveryKthStep == 0) {
			int plot = timeStepsSinceLastPlot / plotEveryKthStep;
			if (plot % 10 == 0)
			{
				printf("%d th plot \n", plot);
				printf("t: %.30f \n", t);
			}
			printParaviewSnapshot(timeStepsSinceLastPlot / plotEveryKthStep);
		}
	}
	end = clock();
	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC; //from https://www.geeksforgeeks.org/how-to-measure-time-taken-by-a-program-in-c/
	printf("Time taken: %lf\n", cpu_time_used);
	closeParaviewVideoFile();

	return 0;
}

//g++ -O3 --std=c++11 summativeCode.c -o summativeCode
//test data:
//2 objects colliding, right heavier, equal velocity:
//summativeCode 10   0 0 0  1 0 0  10   1 0 0  -1 0 0  20
//summativeCode 10.0 0 0 0 1 0 0 10 1 0 0 -1 0 0 20 2 0 0 -0.1 0 0 5
