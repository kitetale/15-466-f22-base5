#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#include <random>

GLuint my_meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer > my_meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("testworld.pnct"));
	my_meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > my_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("testworld.scene"), [&](Scene &scene, Scene::Transform *transform, std::string const &mesh_name){
		Mesh const &mesh = my_meshes->lookup(mesh_name);

		scene.drawables.emplace_back(transform);
		Scene::Drawable &drawable = scene.drawables.back();

		drawable.pipeline = lit_color_texture_program_pipeline;

		drawable.pipeline.vao = my_meshes_for_lit_color_texture_program;
		drawable.pipeline.type = mesh.type;
		drawable.pipeline.start = mesh.start;
		drawable.pipeline.count = mesh.count;

	});
});

WalkMesh const *walkmesh = nullptr;
Load< WalkMeshes > my_walkmeshes(LoadTagDefault, []() -> WalkMeshes const * {
	WalkMeshes *ret = new WalkMeshes(data_path("testworld.w"));
	walkmesh = &ret->lookup("Plane");
	return ret;
});

PlayMode::PlayMode() : scene(*my_scene) {

	//load assests
	//store the leg parts in a vector
	std::vector<Scene::Transform *> leg_parts(6);
	Scene::Transform* body = nullptr;
	for (auto &transform : scene.transforms) {
		if ( transform.name == "leg_l_hip") leg_parts[0] = &transform;
		else if (transform.name == "leg_l_knee") leg_parts[1] = &transform;
		else if (transform.name == "leg_l_ankle") leg_parts[2] = &transform;
		else if (transform.name == "leg_r_hip") leg_parts[3] = &transform;
		else if (transform.name == "leg_r_knee") leg_parts[4] = &transform;
		else if (transform.name == "leg_r_ankle") leg_parts[5] = &transform;
		else if (transform.name == "body") body = &transform;
		else if (transform.name == "Cone") debugcone = &transform;
	}
	if (leg_parts[0] == nullptr) throw std::runtime_error("left hip not found");
	if (leg_parts[1] == nullptr) throw std::runtime_error("left knee not found");
	if (leg_parts[2] == nullptr) throw std::runtime_error("left ankle not found");
	if (leg_parts[3] == nullptr) throw std::runtime_error("right hip not found");
	if (leg_parts[4] == nullptr) throw std::runtime_error("right knee not found");
	if (leg_parts[5] == nullptr) throw std::runtime_error("right ankle not found");
	if (body == nullptr) throw std::runtime_error("body not found");
	if (debugcone == nullptr) throw std::runtime_error("cone not found");
	//construct legs and Walker
	Leg Left = Leg(leg_parts[0], leg_parts[1], leg_parts[2]);
	Leg Right = Leg(leg_parts[3], leg_parts[4], leg_parts[5]);
	std::cout << "does walker fail? \n";
	walker = Walker(body, Left, Right);
	std::cout << "walker made \n";

	//-------------------------------------------------------------
	//create a player transform:
	scene.transforms.emplace_back();
	player.transform = &scene.transforms.back();

	//create a player camera attached to a child of the player transform:
	scene.transforms.emplace_back();
	scene.cameras.emplace_back(&scene.transforms.back());
	player.camera = &scene.cameras.back();
	player.camera->fovy = glm::radians(60.0f);
	player.camera->near = 0.01f;
	player.camera->transform->parent = player.transform;
	//player's eyes are 1.8 units above the ground:
	player.camera->transform->position = glm::vec3(0.0f, 0.0f, 1.8f);

	//rotate camera facing direction (-z) to player facing direction (+y):
	player.camera->transform->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));

	//start player walking at nearest walk point:
	player.at = walkmesh->nearest_walk_point(player.transform->position);

	//initalize Walker's position
	walker.at = walkmesh->nearest_walk_point(walker.body->getWorldPosition());
	walker.body->position = walkmesh->to_world_point(walker.at) + glm::vec3(0.0f, 0.0f, walker.groundOffset);
	//initialize the walkers leg's mesh 
	walker.left_leg.at = walkmesh->nearest_walk_point(walker.left_leg.hip->getWorldPosition() + glm::vec3(1.0f, 0.0f, 0.0f));
	walker.right_leg.at = walkmesh->nearest_walk_point(walker.right_leg.hip->getWorldPosition() + glm::vec3(1.0f, 0.0f, 0.0f));

	//TODO: change into average of normal later, when u set up the leg walk animation
	//snake.left_leg.hip->position = walkmesh->to_world_point(snake.at) + glm::vec3(0.0f, 0.0f, snake.groundOffset);
	//leg_at = walkmesh->nearest_walk_point(right_leg.hip->getWorldPosition());

}

PlayMode::~PlayMode() {
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_ESCAPE) {
			SDL_SetRelativeMouseMode(SDL_FALSE);
			return true;
		} else if (evt.key.keysym.sym == SDLK_a) {
			left.downs += 1;
			left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.downs += 1;
			right.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.downs += 1;
			up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.downs += 1;
			down.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_UP) { //additional debug inputs
			front.downs += 1;
			front.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_DOWN) {
			back.downs += 1;
			back.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_LEFT) {
			turn_left.downs += 1;
			turn_left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_RIGHT) {
			turn_right.downs += 1;
			turn_right.pressed = true;
			return true;
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			up.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			down.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_UP) { //additional debug inputs
			front.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_DOWN) {
			back.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_LEFT) {
			turn_left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_RIGHT) {
			turn_right.pressed = false;
			return true;
		}
	} else if (evt.type == SDL_MOUSEBUTTONDOWN) {
		if (SDL_GetRelativeMouseMode() == SDL_FALSE) {
			SDL_SetRelativeMouseMode(SDL_TRUE);
			return true;
		}
	} else if (evt.type == SDL_MOUSEMOTION) {
		if (SDL_GetRelativeMouseMode() == SDL_TRUE) {
			glm::vec2 motion = glm::vec2(
				evt.motion.xrel / float(window_size.y),
				-evt.motion.yrel / float(window_size.y)
			);
			glm::vec3 upDir = walkmesh->to_world_smooth_normal(player.at);
			player.transform->rotation = glm::angleAxis(-motion.x * player.camera->fovy, upDir) * player.transform->rotation;

			float pitch = glm::pitch(player.camera->transform->rotation);
			pitch += motion.y * player.camera->fovy;
			//camera looks down -z (basically at the player's feet) when pitch is at zero.
			pitch = std::min(pitch, 0.95f * 3.1415926f);
			pitch = std::max(pitch, 0.05f * 3.1415926f);
			player.camera->transform->rotation = glm::angleAxis(pitch, glm::vec3(1.0f, 0.0f, 0.0f));

			return true;
		}
	}

	return false;
}

void PlayMode::update(float elapsed) {

	//testing if leg ik works
	timer += elapsed;
	if (timer > 1.0f) {
		timer = 0.0f;
	}
	
	
	float percent = timer / 5.0f;
	glm::vec3 target = glm::mix(glm::vec3(0.f, 0.0f, 3.0f), glm::vec3(0.f, 0.0f, 0.0f), percent);
	
	
	//player walking:
	{		
		//combine inputs into a move:
		constexpr float PlayerSpeed = 5.0f;
		glm::vec2 move = glm::vec2(0.0f);
		if (left.pressed && !right.pressed) move.x =-1.0f;
		if (!left.pressed && right.pressed) move.x = 1.0f;
		if (down.pressed && !up.pressed) move.y =-1.0f;
		if (!down.pressed && up.pressed) move.y = 1.0f;

		//make it so that moving diagonally doesn't go faster:
		if (move != glm::vec2(0.0f)) move = glm::normalize(move) * PlayerSpeed * elapsed;

		//get move in world coordinate system:
		glm::vec3 remain = player.transform->make_local_to_world() * glm::vec4(move.x, move.y, 0.0f, 0.0f);

		//using a for() instead of a while() here so that if walkpoint gets stuck in
		// some awkward case, code will not infinite loop:
		for (uint32_t iter = 0; iter < 10; ++iter) {
			if (remain == glm::vec3(0.0f)) break;
			WalkPoint end;
			float time;
			walkmesh->walk_in_triangle(player.at, remain, &end, &time);
			player.at = end;
			if (time == 1.0f) {
				//finished within triangle:
				remain = glm::vec3(0.0f);
				break;
			}
			//some step remains:
			remain *= (1.0f - time);
			//try to step over edge:
			glm::quat rotation;
			if (walkmesh->cross_edge(player.at, &end, &rotation)) {
				//stepped to a new triangle:
				player.at = end;
				//rotate step to follow surface:
				remain = rotation * remain;
			} else {
				//ran into a wall, bounce / slide along it:
				glm::vec3 const &a = walkmesh->vertices[player.at.indices.x];
				glm::vec3 const &b = walkmesh->vertices[player.at.indices.y];
				glm::vec3 const &c = walkmesh->vertices[player.at.indices.z];
				glm::vec3 along = glm::normalize(b-a);
				glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
				glm::vec3 in = glm::cross(normal, along);

				//check how much 'remain' is pointing out of the triangle:
				float d = glm::dot(remain, in);
				if (d < 0.0f) {
					//bounce off of the wall:
					remain += (-1.25f * d) * in;
				} else {
					//if it's just pointing along the edge, bend slightly away from wall:
					remain += 0.01f * d * in;
				}
			}
		}

		if (remain != glm::vec3(0.0f)) {
			std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
		}

		//update player's position to respect walking:
		player.transform->position = walkmesh->to_world_point(player.at);

		{ //update player's rotation to respect local (smooth) up-vector:
			
			glm::quat adjust = glm::rotation(
				player.transform->rotation * glm::vec3(0.0f, 0.0f, 1.0f), //current up vector
				walkmesh->to_world_smooth_normal(player.at) //smoothed up vector at walk location
			);
			player.transform->rotation = glm::normalize(adjust * player.transform->rotation);
		}

		/*
		glm::mat4x3 frame = camera->transform->make_local_to_parent();
		glm::vec3 right = frame[0];
		//glm::vec3 up = frame[1];
		glm::vec3 forward = -frame[2];

		camera->transform->position += move.x * right + move.y * forward;
		*/
	}

	static bool once = false;
	//make a lambda to move leg...
	{
		float angle = 0.f;
		float dir = 0.f;
		if (front.pressed && !back.pressed) dir = 1.0f;
		if (back.pressed && !front.pressed) dir = -1.0f;
		if (turn_left.pressed && !turn_right.pressed) angle = 5.0f;
		if (turn_right.pressed && !turn_left.pressed) angle = -5.0f;


		//UPDATING WALKER------------------
		walker.world_rotation += angle; //use this to keep track of leg rotation
		walker.world_rotation = glm::mod(walker.world_rotation, 360.f); //stays within 360
		//apply the world rotation 
		walker.body->rotation = glm::angleAxis(glm::radians(walker.world_rotation), glm::vec3(0.0f, 0.0f, 1.0f));
		//move in facing direction (in world space):
		glm::vec3 offset = glm::rotate(walker.body->rotation, glm::vec3(0.0f, dir*0.05f, 0.0f));
		walker.update_at(offset); //also sets the position
		//now take the legs with it
		walker.update_legs();

		//-------------------------------


		//update Leg
		Leg* leg = &walker.left_leg;
		
			//called in update to do everything necessary to move the leg
		leg->leg_routine(offset, elapsed);
		

		leg = &walker.right_leg;
		leg->update_leg_at(offset);		//update leg position:
		leg->update(leg->get_animated_position());
		//update the step position every time timer restes
		if (leg->is_too_far() && !leg->animating) {
			leg->prev_step = leg->step_to;
			leg->update_step_to();
			leg->animating = true;
		}
		if (leg->animating) {
			leg->timer += elapsed;
			if (leg->timer > leg->anim_time) {
				leg->animating = false;
				leg->timer = 0.f;
			}
		}
		//just rotate 180 around the world z axis
		//walker.body->rotation = glm::angleAxis(glm::radians(180.f), glm::vec3(0.0f, 0.0f, 1.0f)) * walker.body->rotation;
	}
	

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	up.downs = 0;
	down.downs = 0;
	front.downs = 0;
	back.downs = 0;
	turn_left.downs = 0;
	turn_right.downs = 0;
}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	player.camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	//set up light type and position for lit_color_texture_program:
	// TODO: consider using the Light(s) in the scene to do this
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.

	scene.draw(*player.camera);

	/* In case you are wondering if your walkmesh is lining up with your scene, try:
	{
		glDisable(GL_DEPTH_TEST);
		DrawLines lines(player.camera->make_projection() * glm::mat4(player.camera->transform->make_world_to_local()));
		for (auto const &tri : walkmesh->triangles) {
			lines.draw(walkmesh->vertices[tri.x], walkmesh->vertices[tri.y], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
			lines.draw(walkmesh->vertices[tri.y], walkmesh->vertices[tri.z], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
			lines.draw(walkmesh->vertices[tri.z], walkmesh->vertices[tri.x], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
		}
	}
	*/

	{ //use DrawLines to overlay some text:
		glDisable(GL_DEPTH_TEST);
		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		constexpr float H = 0.09f;
		lines.draw_text("Mouse motion looks; WASD moves; escape ungrabs mouse",
			glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		float ofs = 2.0f / drawable_size.y;
		lines.draw_text("Mouse motion looks; WASD moves; escape ungrabs mouse",
			glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0xff, 0xff, 0xff, 0x00));
	}
	GL_ERRORS();
}

void Walker::update_legs() {
	//rotate the leg_to_leg such that it has body's rotation:
	glm::mat4 rot_mat = glm::rotate(glm::mat4(1.0f), glm::radians(360.f - world_rotation), glm::vec3(0.0f, 0.0f, 1.0f));
	glm::vec3 left_offset = glm::vec3(glm::vec4(body_to_leftleg, 1.0f) * rot_mat);
	glm::vec3 right_offset = glm::vec3(glm::vec4(body_to_rightleg, 1.0f) * rot_mat);
	//add onto postion of legs:
	left_leg.hip->position = body->position + left_offset;
	right_leg.hip->position = body->position + right_offset;
	//update the leg's at position 
}
//updates at, and updates position
void Walker::update_at(glm::vec3 remain) {
	for (uint32_t iter = 0; iter < 10; ++iter) {
		if (remain == glm::vec3(0.0f)) break;
		WalkPoint end;
		float time;
		walkmesh->walk_in_triangle(at, remain, &end, &time);
		at = end;
		if (time == 1.0f) {
			//finished within triangle:
			remain = glm::vec3(0.0f);
			break;
		}
		//some step remains:
		remain *= (1.0f - time);
		//try to step over edge:
		glm::quat rotation;
		if (walkmesh->cross_edge(at, &end, &rotation)) {
			//stepped to a new triangle:
			at = end;
			//rotate step to follow surface:
			remain = rotation * remain;
		} else {
			//ran into a wall, bounce / slide along it:
			glm::vec3 const &a = walkmesh->vertices[at.indices.x];
			glm::vec3 const &b = walkmesh->vertices[at.indices.y];
			glm::vec3 const &c = walkmesh->vertices[at.indices.z];
			glm::vec3 along = glm::normalize(b-a);
			glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
			glm::vec3 in = glm::cross(normal, along);

			//check how much 'remain' is pointing out of the triangle:
			float d = glm::dot(remain, in);
			if (d < 0.0f) {
				//bounce off of the wall:
				remain += (-1.25f * d) * in;
			} else {
				//if it's just pointing along the edge, bend slightly away from wall:
				remain += 0.01f * d * in;
			}
		}
	}

	if (remain != glm::vec3(0.0f)) {
		std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
		assert(false);
	}

	//set world position to be at the walk point:
	body->position = walkmesh->to_world_point(at) + glm::vec3(0.0f, 0.0f, groundOffset);

}

//take target (in world space) and updates the ankle to target
void Leg::update(glm::vec3 const &targetWorld) {
	/*referenced this tutorial:
	https://www.alanzucconi.com/2018/05/02/ik-2d-1/ */
	//a little too tired to implement this rn, but here's my thougts:
	
	std::cout << "targetworld: " << targetWorld.x << ", " << targetWorld.y << ", " << targetWorld.z << std::endl;
	glm::vec3 hip_world = hip->position;
	std::cout << "hip_world: " << hip_world.x << ", " << hip_world.y << ", " << hip_world.z << std::endl;
	glm::vec3 omega = targetWorld - hip_world;

	if (abs(omega.x) < 0.0005f) omega.x = 0.f;
	if (abs(omega.y) < 0.0005f) omega.y = 0.0f;
	float theta_base = atan2(omega.y, omega.x);
	glm::quat rotateTo = glm::angleAxis(-theta_base, glm::vec3(0.0f, 0.0f, 1.0f));
	glm::quat rotateBack = glm::angleAxis(theta_base, glm::vec3(0.0f, 0.0f, 1.0f));
	//rotate target such that target is on the x axis
	glm::vec3 target = rotateTo* omega;
	std::cout << "target: " << target.x << ", " << target.y << ", " << target.z << std::endl;
	//glm::quat tester = glm::inverse(body->getWorldRotation()) * rotateBack;

	//rotate the hip back such that it's x axis is the same as the world x axis
	length_c = glm::length(target);// this calculation needs to be done on that plane
	//std::cout << "length_c: " << length_c << std::endl;
	//account for the case where target is too far away:
	//angle from hip to target (in radians)
	float theta = atan2(target.z, target.x);
	//std::cout << "theta: " << theta << std::endl;
	//too far away
	if ( length_a + length_b < length_c) {
		angleA = theta;
		angleB = 0.0f;
	} 
	else
	{
		//float thing = (length_a * length_a + length_c * length_c - length_b * length_b) / (2 * length_a * length_c);
		//std::cout<< "thing: " << thing << std::endl;
		//inner angle Alpha (in radians)
		float alpha = acos((length_a * length_a + length_c * length_c - length_b * length_b) / (2.0f * length_a * length_c));
		//inner angle beta (in radians)
		float beta = acos((length_a * length_a + length_b * length_b - length_c * length_c) / (2.0f * length_a * length_b));
		//calculate the respective angles to be applied to the joints:
		//reversing the signs here can change the direction of the leg
		angleA = theta - alpha;
		angleB = (float)M_PI - beta;
	}	

	//update the transforms:
	hip->rotation = rotateBack * glm::angleAxis(angleA, glm::vec3(0.0f, -1.0f,.0f));
	knee->rotation = glm::angleAxis(angleB, glm::vec3(0.0f, -1.0f, 0.0f));

}

void Leg::printEverything() {
	std::cout << "---------------------" << std::endl;
	//find world position 
	glm::vec3 hipworld = hip->make_local_to_world() * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	glm::vec3 kneeworld = knee->make_local_to_world() * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	glm::vec3 footworld = ankle->make_local_to_world() * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	//print world position
	std::cout << "hipworld: " << hipworld.x << ", " << hipworld.y << ", " << hipworld.z << std::endl;
	std::cout << "kneeworld: " << kneeworld.x << ", " << kneeworld.y << ", " << kneeworld.z << std::endl;
	std::cout << "footworld: " << footworld.x << ", " << footworld.y << ", " << footworld.z << std::endl;
	//so what the fuck are these
	/*
	std::cout << "hip position: " << hip->position.x << ", " << hip->position.y << ", " << hip->position.z << std::endl;
	std::cout << "knee position: " << knee->position.x << ", " << knee->position.y << ", " << knee->position.z << std::endl;
	std::cout << "ankle position: " << ankle->position.x << ", " << ankle->position.y << ", " << ankle->position.z << std::endl;
	*/
	std::cout << "angleA: " << angleA << std::endl;
	std::cout << "angleB: " << angleB << std::endl;
	std::cout << "length_a: " << length_a << std::endl;
	std::cout << "length_b: " << length_b << std::endl;
	std::cout << "---------------------" << std::endl;
}

void Leg::update_leg_at(glm::vec3 remain)
{
	Leg* leg = this;
	for (uint32_t iter = 0; iter < 10; ++iter) {
		if (remain == glm::vec3(0.0f)) break;
		WalkPoint end;
		float time;
		walkmesh->walk_in_triangle(leg->at, remain, &end, &time);
		leg->at = end;
		if (time == 1.0f) {
			//finished within triangle:
			remain = glm::vec3(0.0f);
			break;
		}
		//some step remains:
		remain *= (1.0f - time);
		//try to step over edge:
		glm::quat rotation;
		if (walkmesh->cross_edge(leg->at, &end, &rotation)) {
			//stepped to a new triangle:
			leg->at = end;
			//rotate step to follow surface:
			remain = rotation * remain;
		} else {
			//ran into a wall, bounce / slide along it:
			glm::vec3 const &a = walkmesh->vertices[leg->at.indices.x];
			glm::vec3 const &b = walkmesh->vertices[leg->at.indices.y];
			glm::vec3 const &c = walkmesh->vertices[leg->at.indices.z];
			glm::vec3 along = glm::normalize(b-a);
			glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
			glm::vec3 in = glm::cross(normal, along);

			//check how much 'remain' is pointing out of the triangle:
			float d = glm::dot(remain, in);
			if (d < 0.0f) {
				//bounce off of the wall:
				remain += (-1.25f * d) * in;
			} else {
				//if it's just pointing along the edge, bend slightly away from wall:
				remain += 0.01f * d * in;
			}
		}
	}

	if (remain != glm::vec3(0.0f)) {
		std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
		assert(false);
	}

}

void Leg::update_step_to() {
	step_to = walkmesh->to_world_point(at);
}

bool Leg::is_too_far()
{
	const glm::vec3 removeZ = glm::vec3(1.0f, 1.0f, 0.0f);
	if (glm::distance(step_to * removeZ, hip->position * removeZ) > 5.f) {
		return true;
	}
	return false;
}

glm::vec3 Leg::get_animated_position(){
	if (!animating) return step_to;
	float t = timer/(anim_time);
	glm::vec3 pos = glm::mix(prev_step, step_to, t);
	//make a z offset that goes up and down 
	float zpos = sinf((float)M_PI * t)* 2.0f; //multiply by 2 arbitrary height 
	return pos + glm::vec3(0.0f, 0.0f, zpos);
}

void Leg::leg_routine(glm::vec3 offset, float elapsed)
{
	update_leg_at(offset);		//update leg position:
	update(get_animated_position());
	//update the step position every time timer restes
	if (is_too_far() && !animating) {
		prev_step = step_to;
		update_step_to();
		animating = true;
	}
	if (animating) {
		timer += elapsed;
		if (timer > anim_time) {
			animating = false;
			timer = 0.f;
		}
	}
}