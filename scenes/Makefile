.PHONY : all

#n.b. the '-y' sets autoexec scripts to 'on' so that driver expressions will work
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	BLENDER = /Applications/blender.app/Contents/MacOS/blender -y
else ifeq ($(UNAME_S),Linux)
	BLENDER = ../../blender-2.93.4-linux-x64/blender -y
else
	BLENDER = blender
endif

EXPORT_MESHES=export-meshes.py
EXPORT_WALKMESHES=export-walkmeshes.py
EXPORT_SCENE=export-scene.py

DIST=../dist

all : \
	$(DIST)/maze.pnct \
	$(DIST)/maze.w \
	$(DIST)/maze.scene \

$(DIST)/maze.pnct : maze.blend $(EXPORT_MESHES)
	$(BLENDER) --background --python $(EXPORT_MESHES) -- '$<':Platforms '$@'

$(DIST)/maze.scene : maze.blend $(EXPORT_SCENE)
	$(BLENDER) --background --python $(EXPORT_SCENE) -- '$<':Platforms '$@'

$(DIST)/maze.w : maze.blend $(EXPORT_WALKMESHES)
	$(BLENDER) --background --python $(EXPORT_WALKMESHES) -- '$<':WalkMeshes '$@'
