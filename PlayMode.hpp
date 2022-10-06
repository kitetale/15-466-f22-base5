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

	//for debug purposes, but annoying
	glm::quat hip_start_pos = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
	//walkmesh location:
	WalkPoint at;
	//where it is stepping
	glm::vec3 step_to = glm::vec3(0.0f, 0.0f, 0.0f);
	//prev step(animation)
	glm::vec3 prev_step = glm::vec3(0.0f, 0.0f, 0.0f);
	//timer and animations 
	float timer = 0.0f;
	float anim_time = 0.25f;
	bool animating = false;
	//make a default constructor:
	Leg() = default;
	//make a constructor that takes in the three joints:
	Leg(Scene::Transform *hip_, Scene::Transform *knee_, Scene::Transform *ankle_) : hip(hip_), knee(knee_), ankle(ankle_) {
		//compute lengths of leg segments:
		length_a = glm::length(hip->getWorldPosition() - knee->getWorldPosition());
		length_b = glm::length(knee->getWorldPosition() - ankle->getWorldPosition());
		hip_start_pos = hip->rotation;
	}
	//update the walkmesh position based on current 
	//update function to move the leg to a given position:
	void update(glm::vec3 const &targetWorld);
	//debug print function:
	void printEverything();
	//update leg walk point
	void update_leg_at(glm::vec3 remain);
	//update the stepping point to at position...
	void update_step_to();
	//check if legs are going too far
	bool is_too_far();
	//called in update to do everything necessary to move the leg
	void leg_routine(glm::vec3 offset, float elapsed);
	//returns animated position 
	glm::vec3 get_animated_position();

};


struct Walker {
	//a body with two legs
	Scene::Transform *body = nullptr;
	Leg left_leg = Leg();
	Leg right_leg = Leg();

	//position on the walkmesh:
	WalkPoint at;
	//how far should we hover off the ground?
	float groundOffset = 10.0f;
	//body to leg offset
	glm::vec3 body_to_leftleg = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 body_to_rightleg = glm::vec3(0.0f, 0.0f, 0.0f);

	//movement:
	float speed = 0.25f;

	float world_rotation = 0.0f;

	//make a default constructor:
	Walker() = default;

	//make a constructor that takes in the body and the two legs:
	Walker(Scene::Transform *body_, Leg left_leg_, Leg right_leg_) : body(body_), left_leg(left_leg_), right_leg(right_leg_) {
		body_to_leftleg = left_leg.hip->getWorldPosition() - body->getWorldPosition();
		body_to_rightleg =  right_leg.hip->getWorldPosition() - body->getWorldPosition();
		left_leg.knee->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
		
	}

	//update function to move legs accordingly 
	void update_legs();
	//update leg walk point
	void update_at(glm::vec3 remain);

};

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, down, up, front, back, turn_left, turn_right;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	//player info:
	struct Player {
		WalkPoint at;
		//transform is at player's feet and will be yawed by mouse left/right motion:
		Scene::Transform *transform = nullptr;
		//camera is at player's head and will be pitched by mouse up/down motion:
		Scene::Camera *camera = nullptr;
	} player;

	Walker walker;
	Scene::Transform *debugcone = nullptr;

	float timer = 0.0f;
};
