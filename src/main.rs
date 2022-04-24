use std::{
    sync::mpsc::channel,
    time::{Duration, Instant},
};

use glutin::{
    dpi::LogicalSize,
    event::{Event, VirtualKeyCode, WindowEvent},
    event_loop::ControlFlow,
    Api, GlRequest,
};
use luminance::{
    pipeline::{PipelineState, Viewport},
    render_state::RenderState,
    shader::Uniform,
    UniformInterface, Vertex,
};
use luminance_derive::Semantics;
use luminance_front::{context::GraphicsContext, shader::Program};
use luminance_glutin::GlutinSurface;
use notify::{watcher, DebouncedEvent, RecursiveMode, Watcher};

const SDF_VS: &'static str = include_str!("./shader/sdf.vs");
const SDF_FS: &'static str = include_str!("./shader/sdf.fs");

#[derive(Clone, Copy, Debug, Eq, PartialEq, Semantics)]
pub enum Semantics {
    #[sem(name = "position", repr = "[f32; 2]", wrapper = "VertexPosition")]
    Position,

    #[sem(name = "texcoord", repr = "[f32; 2]", wrapper = "VertexTexcoord")]
    TexCoord,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Vertex)]
#[vertex(sem = "Semantics")]
struct Vertex {
    pos: VertexPosition,
    texcoord: VertexTexcoord,
}

#[derive(Debug, UniformInterface)]
struct ShaderInterface {
    time: Uniform<f32>,
}

const QUAD_VERTICES: [Vertex; 4] = [
    Vertex {
        pos: VertexPosition::new([-1.0, -1.0]),
        texcoord: VertexTexcoord::new([0.0, 1.0]),
    },
    Vertex {
        pos: VertexPosition::new([1.0, -1.0]),
        texcoord: VertexTexcoord::new([1.0, 1.0]),
    },
    Vertex {
        pos: VertexPosition::new([-1.0, 1.0]),
        texcoord: VertexTexcoord::new([0.0, 0.0]),
    },
    Vertex {
        pos: VertexPosition::new([1.0, 1.0]),
        texcoord: VertexTexcoord::new([1.0, 0.0]),
    },
];

const QUAD_INDICES: [u16; 4] = [0, 1, 2, 3];

fn compile_shader(
    surface: &mut GlutinSurface,
    vs_source: &str,
    fs_source: &str,
) -> anyhow::Result<Program<Semantics, (), ShaderInterface>> {
    Ok(surface
        .new_shader_program::<Semantics, (), ShaderInterface>()
        .from_strings(vs_source, None, None, fs_source)?
        .ignore_warnings())
}

fn main() -> anyhow::Result<()> {
    let mut width = 1920u32;
    let mut height = 1080u32;

    let (tx, rx) = channel();

    let mut watcher = watcher(tx, Duration::from_millis(200))?;
    watcher.watch("src/shader", RecursiveMode::Recursive)?;

    let (mut surface, event_loop) = GlutinSurface::new_gl33_from_builders(
        |_, window_builder| {
            window_builder
                .with_inner_size(LogicalSize::new(width, height))
                .with_resizable(true)
        },
        |_, ctx_builder| {
            ctx_builder
                .with_gl(GlRequest::Specific(Api::OpenGl, (3, 3)))
                .with_gl_profile(glutin::GlProfile::Core)
                .with_double_buffer(Some(true))
                .with_vsync(true)
                .with_multisampling(2)
        },
    )?;

    let mut vs_source = String::from(SDF_VS);
    let mut fs_source = String::from(SDF_FS);

    let mut program = compile_shader(&mut surface, &vs_source, &fs_source)?;

    let quad = surface
        .new_tess()
        .set_vertices(&QUAD_VERTICES[..])
        .set_indices(&QUAD_INDICES[..])
        .set_mode(luminance::tess::Mode::TriangleStrip)
        .build()?;

    let mut start = Instant::now();

    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Poll;

        match event {
            Event::LoopDestroyed => {
                return;
            }
            Event::WindowEvent { event, .. } => match event {
                WindowEvent::CloseRequested => {
                    *control_flow = ControlFlow::Exit;
                    return;
                }
                WindowEvent::Resized(size) => {
                    width = size.width;
                    height = size.height;
                }
                WindowEvent::KeyboardInput { input, .. } => match input.virtual_keycode {
                    Some(VirtualKeyCode::R) => {
                        start = Instant::now();
                    }
                    _ => {}
                },
                _ => {}
            },
            Event::MainEventsCleared => {
                match rx.try_recv() {
                    Ok(DebouncedEvent::Write(path)) => {
                        if path.ends_with("sdf.vs") {
                            vs_source = std::fs::read_to_string(&path).unwrap();
                        } else if path.ends_with("sdf.fs") {
                            fs_source = std::fs::read_to_string(&path).unwrap();
                        }

                        match compile_shader(&mut surface, &vs_source, &fs_source) {
                            Ok(new_program) => {
                                println!("Reloaded {}!", path.display());
                                program = new_program;
                            }
                            Err(e) => {
                                eprintln!("{}", e);
                            }
                        }
                    }
                    _ => {}
                }

                let program = &mut program;
                let back_buffer = surface.back_buffer().unwrap();

                surface
                    .new_pipeline_gate()
                    .pipeline(
                        &back_buffer,
                        &PipelineState::default().set_viewport(Viewport::Specific {
                            width,
                            height,
                            x: 0,
                            y: 0,
                        }),
                        |_, mut shd_gate| {
                            shd_gate.shade(program, |mut iface, uni, mut rdr_gate| {
                                iface.set(&uni.time, start.elapsed().as_secs_f32());

                                rdr_gate.render(&RenderState::default(), |mut tess_gate| {
                                    tess_gate.render(&quad)
                                })
                            })
                        },
                    )
                    .assume();

                surface.swap_buffers();
            }
            Event::RedrawRequested(_) => {}
            _ => {}
        }
    });
}
