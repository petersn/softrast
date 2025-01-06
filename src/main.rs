use minifb::{Key, Window, WindowOptions};
use nalgebra::{Matrix4, Point3, Point4, Rotation3, UnitQuaternion, Vector3};

const WIDTH:  usize = 640;
const HEIGHT: usize = 480;

fn model_matrix(
  translation: Vector3<f32>,
  rotation: UnitQuaternion<f32>,
  scale: f32,
) -> Matrix4<f32> {
  let mut model_matrix = Matrix4::identity();
  let rotation = rotation.to_rotation_matrix();
  model_matrix.fixed_view_mut::<3, 3>(0, 0).copy_from(&rotation.matrix());
  model_matrix.fixed_view_mut::<3, 1>(0, 3).copy_from(&translation);
  model_matrix.m43 *= scale;
  model_matrix
}

fn view_matrix(eye: Point3<f32>, look_at: Point3<f32>) -> Matrix4<f32> {
  let view = (look_at - eye).normalize();
  let right = Vector3::y().cross(&view).normalize();
  let up = view.cross(&right);
  Matrix4::new(
    right.x, right.y, right.z, -right.dot(&eye.coords),
    up.x,    up.y,    up.z,    -up.dot(&eye.coords),
    view.x,  view.y,  view.z,  -view.dot(&eye.coords),
    0.0,     0.0,     0.0,      1.0,
  )
}

fn perspective_matrix(
  fov_deg: f32, aspect: f32, near: f32, far: f32,
) -> Matrix4<f32> {
  let fov_rad = fov_deg * std::f32::consts::PI / 180.0;
  let height = 1.0 / (fov_rad / 2.0).tan();
  let width = height * aspect;
  let m22 = far / (far - near);
  Matrix4::new(
    width, 0.0,    0.0,  0.0,
    0.0,   height, 0.0,  0.0,
    0.0,   0.0,    m22, -near * m22,
    0.0,   0.0,    1.0,  0.0,
  )
}

fn project_to_ndc(pvm: &Matrix4<f32>, p: Point3<f32>) -> Point3<f32> {
  let p = pvm * Point4::new(p.x, p.y, p.z, 1.0);
  Point3::new(p.x / p.w, p.y / p.w, p.z / p.w)
}

struct Triangle {
  coords: [Point3<f32>; 3],
  color: u32,
}

fn draw(frame_counter: u32, buf: &mut [u32]) {
  let triangles = [
    Triangle {
      coords: [
        Point3::new(-1.0, -1.0, 0.0),
        Point3::new( 1.0, -1.0, 0.0),
        Point3::new( 0.0,  1.0, 0.0),
      ],
      color: 0x00FF00FF,
    },
  ];
  let t = frame_counter as f32 * 0.01;
  let model_matrix = model_matrix(
    Vector3::new(0.0, 0.0, 0.0),
    Rotation3::from_euler_angles(
      t, t, t,
    ).into(),
    1.0,
  );
  let view_matrix = view_matrix(
    Point3::new(0.0, 0.0, -5.0),
    Point3::new(0.0, 0.0, 0.0),
  );
  let projection_matrix = perspective_matrix(
    60.0, HEIGHT as f32 / WIDTH as f32, 0.1, 100.0,
  );
  let pvm = projection_matrix * view_matrix * model_matrix;

  // Fill with black.
  buf.fill(0x00000000);
  // Draw triangles.
  for triangle in &triangles {
    let ndc = [
      project_to_ndc(&pvm, triangle.coords[0]),
      project_to_ndc(&pvm, triangle.coords[1]),
      project_to_ndc(&pvm, triangle.coords[2]),
    ];
    // Draw a triangle in screen coords.
    for pair in [
      (&ndc[0], &ndc[1]),
      (&ndc[1], &ndc[2]),
      (&ndc[0], &ndc[2]),
    ] {
      for step in 0..200 {
        let t = step as f32 / 200.0;
        let x = (pair.0.x * (1.0 - t) + pair.1.x * t) * WIDTH as f32 / 2.0 + WIDTH as f32 / 2.0;
        let y = (pair.0.y * (1.0 - t) + pair.1.y * t) * HEIGHT as f32 / 2.0 + HEIGHT as f32 / 2.0;
        if x >= 0.0 && x < WIDTH as f32 && y >= 0.0 && y < HEIGHT as f32 {
          let x = x as usize;
          let y = y as usize;
          buf[y * WIDTH + x] = triangle.color;
        }
      }
    }
  }
}

fn main() {
  let mut window = Window::new(
    "Rasterize",
    WIDTH,
    HEIGHT,
    WindowOptions::default(),
  ).unwrap();
  window.set_target_fps(60);

  let mut frame_counter = 0;
  let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];
  while window.is_open() && !window.is_key_down(Key::Escape) {
    draw(frame_counter, &mut buffer);
    window.update_with_buffer(&buffer, WIDTH, HEIGHT).unwrap();
    frame_counter += 1;
  }
}
