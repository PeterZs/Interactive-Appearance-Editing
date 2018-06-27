// standard global variables
var scene, camera, renderer, textureLoader, light;

// Character 3d object
var character = null;

// FUNCTIONS
function init() {

  // SCENE
  scene = new THREE.Scene();
  textureLoader = new THREE.TextureLoader();

  // CAMERA
  var SCREEN_WIDTH = window.innerWidth;
  var SCREEN_HEIGHT = window.innerHeight;
  var VIEW_ANGLE = 45;
  var ASPECT = SCREEN_WIDTH / SCREEN_HEIGHT;
  var NEAR = 0.1;
  var FAR = 1000;
  camera = new THREE.PerspectiveCamera(VIEW_ANGLE, ASPECT,
                                       NEAR, FAR);
  scene.add(camera);
  camera.position.set(0,0,5);
  camera.lookAt(scene.position);

  // RENDERER
  renderer = new THREE.WebGLRenderer({
    antialias:true,
    alpha: true
  });
  renderer.setSize(SCREEN_WIDTH, SCREEN_HEIGHT);
  var container = document.body;
  container.appendChild( renderer.domElement );

  // Create light
  light = new THREE.PointLight(0xffffff, 1.0);
  // We want it to be very close to our character
  light.position.set(0.0, 0.0, 0.1);
  scene.add(light);

  // Create character
  character = buildCharacter();
  scene.add(character);

  // Start animation
  animate();
}

var buildCharacter = (function() {
  var _geo = null;

  // Share the same geometry across all planar objects
  function getPlaneGeometry() {
    if(_geo == null) {
      _geo = new THREE.SphereGeometry(1.0, 1.0);
    }

    return _geo;
  };

  return function() {
    var g = getPlaneGeometry();

    var creatureImage = textureLoader.load('assets/textures/basketball.png');
    creatureImage.magFilter = THREE.NearestFilter;

    var mat = new THREE.ShaderMaterial({
        uniforms: THREE.UniformsUtils.merge([
            THREE.UniformsLib['lights'],
            {
                lightIntensity: {type: 'f', value: 1.0},
                textureSampler: {type: 't', value: null}
            }
        ]),
        vertexShader: document.
                      getElementById('vertShader').text,
        fragmentShader: document.
                        getElementById('fragShader').text,
        transparent: true,
        lights: true
    });
    // THREE.UniformsUtils.merge() call THREE.clone() on
    // each uniform. We don't want our texture to be
    // duplicated, so I assign it to the uniform value
    // right here.
    mat.uniforms.textureSampler.value = creatureImage;

    var obj = new THREE.Mesh(g, mat);

    return obj;
  }
})();

function animate() {
  // Update light profile
  var timestampNow = new Date().getTime()/1000.0;
  var lightIntensity = 0.75 +
                       0.25 * Math.cos(timestampNow *
                                       Math.PI);

  character.material.uniforms
           .lightIntensity.value = lightIntensity;
  light.color.setHSL(lightIntensity, 1.0, 0.5);

  // Render scene
  renderer.render(scene, camera);
  requestAnimationFrame(animate);
}