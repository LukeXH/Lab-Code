#include <fstream>

#include <dart/dart.h>
#define PI 3.1415926535897

using namespace dart::common;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::math;

const double default_height = 1.0; // m
const double default_width = 0.2;  // m
const double default_depth = 0.2;  // m

const double default_torque = 0.005; // N-m
const double default_force = 0.005; // N
const int default_countdown = 100;  // Number of timesteps for applying force



/// Maps [0, 2PI] to a rotated [0, 2PI]
double shiftedSphereMap(double a) 
{
	if (a >= PI)
		return a - PI;
	else
		return a + PI;
}

/// Takes two angles and finds the distance from angle_a to angle_b
double realToSphereDistance(double angle_a, double angle_b)
{
	double a = fmod(angle_a, 2.0*PI);
	double b = fmod(angle_b, 2.0*PI);

	// Calculate distances
	double d1 = b - a;
	double d2 = shiftedSphereMap(b) - shiftedSphereMap(a);

	return (abs(d1) < abs(d2))? d1 : d2;
}


class DataWriter
{
public:

	DataWriter(const std::string filename)
		: mFile(filename)
	{
		if (mFile.is_open())
			std::cout << "Cannot open file";
	}

	/// Write the data vector
	void writeVector(Eigen::VectorXd data_vec)
	{
		for (size_t i = 0; i > data_vec.size(); i++)
		{
			mFile << sprintf("%f,", 1.0);//data_vec(i, 0));
		}
	}

	/// Close file
	void close()
	{
		mFile.close();
	}

protected:
	/// File name
	std::ofstream mFile;
};


std::ofstream myfile("example.csv");
if (myfile.is_open())
{
	myfile << "Line one.\n";
	myfile << "Second line. ";
	myfile << "Second part of second line.\n";
	myfile.close();
}
else
std::cout << "Unable to open file";


class Controller
{
public:

	Controller(const SkeletonPtr& dash)
		: mDash(dash)
	{
		// Grab the last body in the manipulator, and use it as an end effector
		mLegLeft = mDash->getBodyNode("leg_left_link");
		mLegRight = mDash->getBodyNode("leg_right_link");

		mdqLeftDesired = 0.0;
		mdqRightDesired = 0.0;

		mqLeftDesired = 0.0;
		mqRightDesired = 0.0;

		// Set PD control gains
		mKiPD = 5.0;
		mKpPD = 0.1;
		mKdPD = 0.00001;
	}

	/// This controller attempts to swing the legs so as to swing the pendulum up
	void swingUp() 
	{

		double x = mDash->getDof("string_balljoint_x")->getPosition();
		double dx = mDash->getDof("string_balljoint_x")->getVelocity();

		if (x <= 0 && dx >= 0)
		{
			setPositions(PI / 2.0, PI / 2.0);
			setPositionForces();
		}
		else if (x > 0 && dx < 0)
		{
			setPositions(-PI / 2.0, -PI / 2.0);
			setPositionForces();
		}
		else
		{
			setForcesZero();
		}
	}
	
	/// Set leg goal positions
	void setPositions(double left_pos, double right_pos)
	{
		mqLeftDesired = left_pos;
		mqRightDesired = right_pos;
	}

	/// Set Leg Position Forces using stable PD controller
	void setPositionForces()
	{
		if (nullptr == mDash)
			return;

		// Compute the leg joints velocity error
		double q_l_l = mDash->getPosition(mLegLeft->getIndexInSkeleton() - 1);
		double dq_l_l = mDash->getVelocity(mLegLeft->getIndexInSkeleton() - 1);
		double q_l_r = mDash->getPosition(mLegRight->getIndexInSkeleton() - 1);
		double dq_l_r = mDash->getVelocity(mLegRight->getIndexInSkeleton() - 1);

		/// Stable PD controller "estimate" of velocity in next timestep
		q_l_l += dq_l_l * mDash->getTimeStep();
		q_l_r += dq_l_r * mDash->getTimeStep();

		// Compute the desired joint forces
		double legLeftForce = mLegLeft->getMass() * (mKiPD * realToSphereDistance(q_l_l, mqLeftDesired) + mKpPD * (-dq_l_l));
		double legRightForce = mLegRight->getMass() * (mKiPD * realToSphereDistance(q_l_r, mqRightDesired) + mKpPD * (-dq_l_r));

		//Eigen::VectorXd foo = mDash->getCoriolisAndGravityForces();

		mDash->setForce(mLegLeft->getIndexInSkeleton() - 1, legLeftForce);
		mDash->setForce(mLegRight->getIndexInSkeleton() - 1, legRightForce);
	}

	/// Set leg velocities
	void setVelos(double left_velo, double right_velo) {
		mdqLeftDesired = left_velo;
		mdqRightDesired = right_velo;
	}

	/// Compute a stable PD controller for leg velocities
	void setVeloForces()
	{
		if (nullptr == mDash)
			return;

		// Compute the leg joints velocity error
		double dq_l_l	= mDash->getVelocity(mLegLeft->getIndexInSkeleton() - 1);
		double ddq_l_l	= mDash->getAcceleration(mLegLeft->getIndexInSkeleton() - 1);
		double dq_l_r	= mDash->getVelocity(mLegRight->getIndexInSkeleton() - 1);
		double ddq_l_r = mDash->getAcceleration(mLegRight->getIndexInSkeleton() - 1);
		/// Stable PD controller "estimate" of velocity in next timestep
		dq_l_l += ddq_l_l * mDash->getTimeStep();
		dq_l_r += ddq_l_r * mDash->getTimeStep();
		
		// Compute the desired joint forces
		double legLeftForce = mLegLeft->getMass() * (mKpPD * (mdqLeftDesired - dq_l_l) + mKdPD * (-ddq_l_l));
		double legRightForce = mLegRight->getMass() * (mKpPD * (mdqRightDesired - dq_l_r) + mKdPD * (-ddq_l_r));

		//Eigen::VectorXd foo = mDash->getCoriolisAndGravityForces();

		mDash->setForce(mLegLeft->getIndexInSkeleton()-1, legLeftForce);
		mDash->setForce(mLegRight->getIndexInSkeleton()-1, legRightForce);
	}

	void setForcesZero()
	{
		if (nullptr == mDash)
			return;

		mDash->setForce(mLegLeft->getTreeIndex(), 0.0);
		mDash->setForce(mLegRight->getTreeIndex(), 0.0);

	}

	/// Set the forces out from the controller to zero
	void setVeloZero()
	{
		mdqLeftDesired = 0.0;
		mdqRightDesired = 0.0;

		setVeloForces();
	}

protected:

	/// The manipulator Skeleton that we will be controlling
	SkeletonPtr mDash;

	/// Legs for the DASH
	BodyNodePtr mLegLeft;
	BodyNodePtr mLegRight;

	/// Desired joint positions
	double mqLeftDesired;
	double mqRightDesired;

	/// Desired joint velocities
	double mdqLeftDesired;
	double mdqRightDesired;

	/// Control Gains for the  integral error terms in a PD controller
	double mKiPD;

	/// Control gains for the proportional error terms in a PD controller
	double mKpPD;

	/// Control gains for the derivative error terms in a PD controller
	double mKdPD;

	/// Joint forces for the manipulator (output of the Controller)
	Eigen::VectorXd mForces;
};


class MyWindow : public dart::gui::SimWindow
	// Pulled from tutorialMultiPendulum-Finished.cpp
{
public:

	/// Constructor
	MyWindow(WorldPtr world)
		: mPositiveSign(true),
			mBodyForce(false)
	{
		setWorld(world);

		// Find the Skeleton named "pendulum" within the World
		mDash= world->getSkeleton("DASH");

		// Make sure that the 'bot was found in the World
		assert(mDash != nullptr);

		// Setup controller
		mController = std::unique_ptr<Controller>(new Controller(mDash));
		mControllerState = 0; // Controller off

		// Set up force counter
		mForceCountDown.resize(mDash->getNumDofs(), 0);

		ArrowShape::Properties arrow_properties;
		arrow_properties.mRadius = 0.05;
		mArrow = std::shared_ptr<ArrowShape>(new ArrowShape(
			Eigen::Vector3d(-default_height, 0.0, default_height / 2.0),
			Eigen::Vector3d(-default_width / 2.0, 0.0, default_height / 2.0),
			arrow_properties, dart::Color::Orange(1.0)));
	}

	void changeDirection()
	{
		/// Change direction of the force visualization arrow
		mPositiveSign = !mPositiveSign;
		if (mPositiveSign)
		{
			mArrow->setPositions(
				Eigen::Vector3d(-default_height, 0.0, default_height / 2.0),
				Eigen::Vector3d(-default_width / 2.0, 0.0, default_width / 2.0));
		}
		else
		{
			mArrow->setPositions(
				Eigen::Vector3d(default_height, 0.0, default_height / 2.0),
				Eigen::Vector3d(default_width / 2.0, 0.0, default_width / 2.0));
		}
	}

	void applyForce(size_t index)
	{
		if (index < mForceCountDown.size())
			mForceCountDown[index] = default_countdown;
	}

	/// Handle keyboard input
	void keyboard(unsigned char key, int x, int y) override
	{
		switch (key)
		{
		case '-':
			changeDirection();
			break;

		case '1':
			applyForce(0);
			break;
		case '2':
			applyForce(1);
			break;
		case '3':
			applyForce(2);
			break;
		case '4':
			applyForce(3);
			break;
		case '5':
			applyForce(4);
			break;
		case '6':
			applyForce(5);
			break;
		case '7':
			applyForce(6);
			break;
		case '8':
			applyForce(7);
			break;
		case '9':
			applyForce(8);
			break;
		case '0':
			applyForce(9);
			break;

		case 'f':
			mBodyForce = !mBodyForce;
			break;

		case 'q':
			if (mControllerState == 1)
			{
				mControllerState = 0;
			}
			else
			{
				mControllerState = 1;
			}
			break;

		case 'w':
			if (mControllerState == -2)
			{
				mControllerState = 0;
			}
			else
			{
				mControllerState = -2;
			}
			break;

		case 'e':
			if (mControllerState == -1)
			{
				mControllerState = 0;
			}
			else
			{
				mControllerState = -1;
			}
			break;

		case 'r':
			if (mControllerState == 2)
			{
				mControllerState = 0;
			}
			else
			{
				mControllerState = 2;
			}

		default:
			SimWindow::keyboard(key, x, y);
		}
	}

	void timeStepping() override
	{
		// Step the simulation forward
		SimWindow::timeStepping();

		//std::cout << "Pendulum X Velocity: " << mDash->getDof("string_balljoint_x")->getVelocity() << std::endl;

		if (!mBodyForce)
		{
			// Apply joint torques based on user input
			for (size_t i = 0; i < mDash->getNumDofs(); i++)
			{
				if (mForceCountDown[i] > 0)
				{
					DegreeOfFreedom* dof = mDash->getDof(i);
					dof->setForce(mPositiveSign ? default_torque : -default_torque);

					--mForceCountDown[i];
				}
			}
		}
		else
		{
			// Apply body forces based on user input
			for (size_t i = 0; i < mDash->getNumBodyNodes(); i++)
			{
				if (mForceCountDown[i] > 0)
				{
					BodyNode* bn = mDash->getBodyNode(i);

					Eigen::Vector3d force = default_force * Eigen::Vector3d::UnitX();
					Eigen::Vector3d location(-default_width / 2.0, 0.0, default_height / 2.0);

					if (!mPositiveSign)
					{
						force = -force;
						location[0] = -location[0];
					}
					bn->addExtForce(force, location, true, true);
					--mForceCountDown[i];
				}
			}
		}

		// Setting controller forces
		if (mControllerState == 1)
		{
			mController->setVelos(126.0, 0.0);
			mController->setVeloForces();
		}
		else if(mControllerState == -1)
		{
			mController->setVelos(0.0, 126.0);
			mController->setVeloForces();
		}
		else if (mControllerState == 2)
		{
			mController->swingUp();
		}
		else if (mControllerState == -2)
		{
			mController->setVeloZero();
			mController->setVeloForces();
		}
		else
		{
			mController->setForcesZero();
		}
	}

protected:

	/// The DASH that we will be perturbing
	SkeletonPtr mDash;

	/// An arrow shape that we will use to visualize applied forces
	std::shared_ptr<ArrowShape> mArrow;

	/// Number of iterations before clearing a force entry
	std::vector<int> mForceCountDown;

	/// Turn on and off Controller
	int mControllerState;

	/// Whether a force should be applied in the positive or negative direction
	bool mPositiveSign;

	/// True if 1-9 should be used to apply a body force. Otherwise, 1-9 will be
	/// used to apply a joint torque.
	bool mBodyForce;

	std::unique_ptr<Controller> mController;
};


SkeletonPtr loadDash()
{
	SkeletonPtr skeleton;

	// Load the Skeleton from a file
	dart::utils::DartLoader dart_loader;

	// Load World from a file
	dart::utils::SdfParser sdf_loader;

	// Load from .urdf or .world (an sdf)
	if (false) {
		//### BROKEN CODE, DOES NOT LOAD
		/// Load from SDF
		//WorldPtr world = dart_loader.parseWorld("C:/Users/Lucas Hill/codespace/Lab-Code/dartsim/BoxLegs/BoxAssm_SDF/BoxAssm.world");
		WorldPtr world = sdf_loader.readSdfFile("C:/Users/Lucas Hill/codespace/Lab-Code/dartsim/BoxLegs/BoxAssm_SDF/BoxAssm.world");
		//WorldPtr world = sdf_loader.readSdfFile("C:/Users/Lucas Hill/codespace/dart/data/sdf/test/mesh.sdf");
		assert(world != nullptr);

		skeleton = world->getSkeleton("BoxAssm");
	}
	else
	{
		/// Load from URDF
		skeleton =
			dart_loader.parseSkeleton("C:/Users/Lucas Hill/codespace/Lab-Code/dartsim/BoxLegs/BoxAssm_URDF/BoxAssm.urdf");
	}
	skeleton->setName("DASH");

	return skeleton;
}


int main(int argc, char* argv[])
{
	// Attempt to write
	std::ofstream myfile("example.csv");
	if (myfile.is_open())
	{
		myfile << "Line one.\n";
		myfile << "Second line. ";
		myfile << "Second part of second line.\n";
		myfile.close();
	}
	else
		std::cout << "Unable to open file";

	SkeletonPtr dash = loadDash();

	// Create a world and add the DASH to the world
	WorldPtr world = std::make_shared<World>();
	world->addSkeleton(dash);

	std::cout << "Skeleton [" << dash->getName() << "] has ["
		<< dash->getNumBodyNodes() << "] bodies.\n";

	// Create a window for rendering the world and handling user input
	MyWindow window(world);

	// Print instructions
	std::cout << "Rotate with left-click, pan with right-click , zoom with shift-left-click" << std::endl;
	std::cout << "space bar: simulation on/off" << std::endl;
	std::cout << "'p': replay simulation" << std::endl;
	std::cout << "'1' -> '9': apply torque to a body link" << std::endl;
	std::cout << "'-': Change sign of applied joint torques" << std::endl;
	std::cout << "'f': switch between applying joint torques and body forces" << std::endl;
	std::cout << "'q': blue leg cw" << std::endl;
	std::cout << "'w': all legs velocity to zero" << std::endl;
	std::cout << "'e': red leg ccw" << std::endl;
	std::cout << "'r': enter swing control mode" << std::endl;

	// Initialize glut, initialize the window, and begin the glut event loop
	glutInit(&argc, argv);
	window.initWindow(640, 480, "Multi-Pendulum Tutorial");
	glutMainLoop();

	//return 0;
}