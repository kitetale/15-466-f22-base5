#include "Mode.hpp"

#include "Scene.hpp"
#include "WalkMesh.hpp"

#include <glm/glm.hpp>

#include <vector>
#include <deque>

//testing IK: leg
struct Leg {
	//basic transformations:
	Scene::Transform *hip = nullptr;
	Scene::Transform *knee = nullptr;
	Scene::Transform *ankle = nullptr;
	//length of leg segments:
	float length_a = 0.0f;
	float length_b = 0.0f;
	float length_c = 0.0f;
	//supplementary angles for IK:
	float angleA = 0.0f;
	float angleB = 0.0f;

	//walkmesh location:
	WalkPoint at;
	//make a default constructor:

	float outer_angle = 0.0f;
	Leg() = default;
	//make a constructor that takes in the three joints:
	Leg(Scene::Transform *hip_, Scene::Transform *knee_, Scene::Transform *ankle_) : hip(hip_), knee(knee_), ankle(ankle_) {
		//compute lengths of leg segments:
		length_a = glm::length(hip->getWorldPosition() - knee->getWorldPosition());
		length_b = glm::length(knee->getWorldPosition() - ankle->getWorldPosition());
	}

	//update function to move the leg to a given position:
	void update(glm::vec3 const &targetWorld);
	//debug print function:
	void printEverything();
};


struct Walker {
	//a body with two legs
	Scene::Transform *body = nullptr;
	Leg left_leg = Leg();
	Leg right_leg = Leg();

	//position on the walkmesh:
	WalkPoint at;
	//how far should we hover off the ground?
	float groundOffset = 8.0f;
	//body to leg offset:
	glm::vec3 body_to_leftleg = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 body_to_rightleg = glm::vec3(0.0f, 0.0f, 0.0f);

	//movement:
	float speed = 0.5f;

	float world_rotation = 0.0f;

	//make a default constructor:
	Walker() = default;

	//make a constructor that takes in the body and the two legs:
	Walker(Scene::Transform *body_, Leg left_leg_, Leg right_leg_) : body(body_), left_leg(left_leg_), right_leg(right_leg_) {
		body_to_leftleg = body->getWorldPosition() - left_leg.hip->getWorldPosition();
		body_to_rightleg = body->getWorldPosition() - right_leg.hip->getWorldPosition();
		left_leg.knee->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
		
	}

	//update function to move legs accordingly 
	void update_legs();
};

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	//----- game state -----
	//collision check:
	struct Ray {
		glm::vec3 orig;
		glm::vec3 dir;
	} r1, r2, r3, r4;

	glm::vec3 r1_base;

	//width, length, height of character
	glm::vec3 dimension = glm::vec3 (1.0f, 1.0f, 2.0f); // height is from the floor
	glm::vec3 halfDim = dimension/2.0f; //half of each for ease of computation
	glm::vec3 minBound = glm::vec3 (0.0f, 0.0f, 0.0f); // in world coord, update by player pos
	glm::vec3 maxBound = glm::vec3 (1.0f, 1.0f, 2.0f); // in world coord, update by player pos

	bool BoxRayCollision(Ray r);

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	Scene::Transform *ray1 = nullptr;
	glm::quat ray1_base_rot;

	float wobble = 0.0f;

	//player info:
	struct Player {
		WalkPoint at;
		//transform is at player's feet and will be yawed by mouse left/right motion:
		Scene::Transform *transform = nullptr;
	} player;

	// camera 
	struct Character {
		Scene::Transform *character_transform = nullptr;
		Scene::Camera *camera = nullptr;
	} character;
	glm::quat char_base_rot;
};
