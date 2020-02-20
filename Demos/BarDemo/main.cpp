#include "Common/Common.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/TimeStepController.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Simulation/Simulation.h"

#include "helpers.hpp"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

void initParameters();
void timeStep();
void buildModel();
void createMesh();
void render();
void reset(TrajectoryData*);

const unsigned int width = 10;
const unsigned int height = 10;
const unsigned int depth = 50;

const Real defaultDimensionScaler = 0.3;
const Real dimensionScaler = 0.3;
const Real massScaler = 0.5;
const Real singleParticleMass = massScaler * (1.0 / defaultDimensionScaler / defaultDimensionScaler / defaultDimensionScaler) * dimensionScaler * dimensionScaler * dimensionScaler;
// const Real singleParticleMass = 1.0;

const unsigned int c_width = width;
const unsigned int c_height = height;
const unsigned int c_depth = 1;

short simulationMethod = 2;

void initPositions(TrajectoryData *trajServer, SimulationModel *model)
{
	// initial all particles to the initial pose
	ParticleData &pd = model->getParticles();

	for (unsigned int i = 0; i < width; i++)
	{
		for (unsigned int j = 0; j < height; j++)
		{
			for (unsigned int k = 0; k < depth; k++)
			{
				unsigned int idx = i * height * depth + j * depth + k;
				auto init_pos = pd.getPosition0(idx);
				Vector3r new_pos = init_pos;
				trajServer->transform(new_pos);
				if (i == 0 && j == 0)
				{
					// std::cout << init_pos << std::endl;
					// std::cout << new_pos << std::endl;
				}
				pd.setPosition(idx, new_pos);
				pd.setVelocity(idx, Vector3r(0, 0, 0));
				pd.setAcceleration(idx, Vector3r(0, 0, 0));
			}
		}
	}
}

bool updatePositions(TrajectoryData *trajServer, SimulationModel *model)
{
	// update static particles according to the trajectory
	ParticleData &pd = model->getParticles();
	auto currentTime = TimeManager::getCurrent()->getTime();

	if (trajServer->needUpdate(currentTime))
	{
		for (unsigned int i = 0; i < c_width; i++)
		{
			for (unsigned int j = 0; j < c_height; j++)
			{
				for (unsigned int k = 0; k < c_depth; k++)
				{
					unsigned int idx = i * height * depth + j * depth + k;
					auto init_pos = pd.getPosition0(idx);
					Vector3r new_pos = init_pos;
					trajServer->transform(new_pos);
					pd.setPosition(idx, new_pos);
				}
			}
		}
		trajServer->step(currentTime);
		return true;
	}
	else
	{
		return false;
	}
}

Vector3r getPosition(SimulationModel* model, TrajectoryData *trajServer)
{
	// get the current position of the target point (0, -h/2, d)
	int i = 0;
	int j = (width - 1) / 2.0;
	int k = depth - 1;
	int targetIdx = i * height * depth + j * depth + k;

	ParticleData &pd = model->getParticles();
	Vector3r targetPosition = pd.getPosition(targetIdx);
	
	// Transform back to base coordinate
	Eigen::Isometry3d lastTransform_inv = trajServer->getLastTransform().inverse();
	targetPosition = lastTransform_inv * targetPosition;

	return targetPosition;
}

void sim_init(int argc, char **argv)
{
	SimulationModel *model = new SimulationModel();
	model->init();
	Simulation::getCurrent()->setModel(model);
	buildModel();
}

Vector3r sim_exec(std::vector<Eigen::Isometry3d> transforms, Eigen::VectorXd dt_inv)
{
	SimulationModel *model = Simulation::getCurrent()->getModel();

	TrajectoryData *trajServer = new TrajectoryData(transforms, dt_inv);
	initPositions(trajServer, model);

	reset(trajServer);

	const unsigned int numStepsPerItr = 1;
	while (!trajServer->finished())
	{
		for (unsigned int i=0; i< numStepsPerItr; i++)
		{
			START_TIMING("SimStep");

			Simulation::getCurrent()->getTimeStep()->step(*model);
			updatePositions(trajServer, model);

			STOP_TIMING_AVG;
		}
		for (unsigned int i = 0; i < model->getTetModels().size(); i++)
		{
			model->getTetModels()[i]->updateMeshNormals(model->getParticles());
		}
	}
	return getPosition(model, trajServer);
}

int main(int argc, char **argv)
{
	sim_init(argc, argv);
	TrajectoryData fakeInput("../Demos/BarDemo/example_traj_data/");
	auto targetPosition = sim_exec(fakeInput.transforms, fakeInput.dt_inv);
	// sim_exec(fakeInput.transforms, fakeInput.dt_inv);
	// sim_exec(fakeInput.transforms, fakeInput.dt_inv);
	// sim_exec(fakeInput.transforms, fakeInput.dt_inv);

	std::cout << "TARGET COORDINATE: " << targetPosition[0] << " " << targetPosition[1] << " " << targetPosition[2] << std::endl;
}

void reset(TrajectoryData *trajServer)
{
	trajServer->reset();
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Simulation::getCurrent()->reset();

	Simulation::getCurrent()->getModel()->cleanup();
	buildModel();
}

void buildModel()
{
	TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));

	createMesh();
}

void createMesh()
{
	Vector3r offset = -1.0 * Vector3r((Real)(width - 1) / 2.0, (Real)(height - 1) / 2.0, (Real)0.0);
	Vector3r points[width * height * depth];
	for (unsigned int i = 0; i < width; i++)
	{
		for (unsigned int j = 0; j < height; j++)
		{
			for (unsigned int k = 0; k < depth; k++)
			{
				auto &targetPnt = points[i * height * depth + j * depth + k];
				targetPnt = Vector3r((Real)i, (Real)j, (Real)k) + offset;
				targetPnt *= dimensionScaler;
			}
		}
	}

	vector<unsigned int> indices;
	for (unsigned int i = 0; i < width - 1; i++)
	{
		for (unsigned int j = 0; j < height - 1; j++)
		{
			for (unsigned int k = 0; k < depth - 1; k++)
			{
				// For each block, the 8 corners are numerated as:
				//     4*-----*7
				//     /|    /|
				//    / |   / |
				//  5*-----*6 |
				//   | 0*--|--*3
				//   | /   | /
				//   |/    |/
				//  1*-----*2
				unsigned int p0 = i * height * depth + j * depth + k;
				unsigned int p1 = p0 + 1;
				unsigned int p3 = (i + 1) * height * depth + j * depth + k;
				unsigned int p2 = p3 + 1;
				unsigned int p7 = (i + 1) * height * depth + (j + 1) * depth + k;
				unsigned int p6 = p7 + 1;
				unsigned int p4 = i * height * depth + (j + 1) * depth + k;
				unsigned int p5 = p4 + 1;

				// Ensure that neighboring tetras are sharing faces
				if ((i + j + k) % 2 == 1)
				{
					indices.push_back(p2);
					indices.push_back(p1);
					indices.push_back(p6);
					indices.push_back(p3);
					indices.push_back(p6);
					indices.push_back(p3);
					indices.push_back(p4);
					indices.push_back(p7);
					indices.push_back(p4);
					indices.push_back(p1);
					indices.push_back(p6);
					indices.push_back(p5);
					indices.push_back(p3);
					indices.push_back(p1);
					indices.push_back(p4);
					indices.push_back(p0);
					indices.push_back(p6);
					indices.push_back(p1);
					indices.push_back(p4);
					indices.push_back(p3);
				}
				else
				{
					indices.push_back(p0);
					indices.push_back(p2);
					indices.push_back(p5);
					indices.push_back(p1);
					indices.push_back(p7);
					indices.push_back(p2);
					indices.push_back(p0);
					indices.push_back(p3);
					indices.push_back(p5);
					indices.push_back(p2);
					indices.push_back(p7);
					indices.push_back(p6);
					indices.push_back(p7);
					indices.push_back(p0);
					indices.push_back(p5);
					indices.push_back(p4);
					indices.push_back(p0);
					indices.push_back(p2);
					indices.push_back(p7);
					indices.push_back(p5);
				}
			}
		}
	}
	SimulationModel *model = Simulation::getCurrent()->getModel();
	model->addTetModel(width * height * depth, (unsigned int)indices.size() / 4u, points, indices.data());

	ParticleData &pd = model->getParticles();
	for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
	{
		pd.setMass(i, singleParticleMass);
	}
	for (unsigned int i = 0; i < c_width; i++)
	{
		for (unsigned int j = 0; j < c_height; j++)
		{
			for (unsigned int k = 0; k < c_depth; k++)
				pd.setMass(i * height * depth + j * depth + k, 0.0);
		}
	}

	// init constraints
	for (unsigned int cm = 0; cm < model->getTetModels().size(); cm++)
	{
		const unsigned int nTets = model->getTetModels()[cm]->getParticleMesh().numTets();
		const unsigned int *tets = model->getTetModels()[cm]->getParticleMesh().getTets().data();
		const IndexedTetMesh::VertexTets *vTets = model->getTetModels()[cm]->getParticleMesh().getVertexTets().data();
		if (simulationMethod == 1)
		{
			const unsigned int offset = model->getTetModels()[cm]->getIndexOffset();
			const unsigned int nEdges = model->getTetModels()[cm]->getParticleMesh().numEdges();
			const IndexedTetMesh::Edge *edges = model->getTetModels()[cm]->getParticleMesh().getEdges().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const unsigned int v1 = edges[i].m_vert[0] + offset;
				const unsigned int v2 = edges[i].m_vert[1] + offset;

				model->addDistanceConstraint(v1, v2);
			}

			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v1 = tets[4 * i];
				const unsigned int v2 = tets[4 * i + 1];
				const unsigned int v3 = tets[4 * i + 2];
				const unsigned int v4 = tets[4 * i + 3];

				model->addVolumeConstraint(v1, v2, v3, v4);
			}
		}
		else if (simulationMethod == 2)
		{
			TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v1 = tets[4 * i];
				const unsigned int v2 = tets[4 * i + 1];
				const unsigned int v3 = tets[4 * i + 2];
				const unsigned int v4 = tets[4 * i + 3];

				model->addFEMTetConstraint(v1, v2, v3, v4);
			}
		}
		else if (simulationMethod == 3)
		{
			TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v1 = tets[4 * i];
				const unsigned int v2 = tets[4 * i + 1];
				const unsigned int v3 = tets[4 * i + 2];
				const unsigned int v4 = tets[4 * i + 3];

				model->addStrainTetConstraint(v1, v2, v3, v4);
			}
		}
		else if (simulationMethod == 4)
		{
			TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
			for (unsigned int i = 0; i < nTets; i++)
			{
				const unsigned int v[4] = {tets[4 * i], tets[4 * i + 1], tets[4 * i + 2], tets[4 * i + 3]};
				// Important: Divide position correction by the number of clusters
				// which contain the vertex.
				const unsigned int nc[4] = {vTets[v[0]].m_numTets, vTets[v[1]].m_numTets, vTets[v[2]].m_numTets, vTets[v[3]].m_numTets};
				model->addShapeMatchingConstraint(4, v, nc);
			}
		}
		model->getTetModels()[cm]->updateMeshNormals(pd);
	}

	LOG_INFO << "Number of tets: " << indices.size() / 4;
	LOG_INFO << "Number of vertices: " << width * height * depth;
	LOG_INFO << "Finished creating the msh";
}
