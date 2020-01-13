use std::{io, fs, cmp, thread};
use image::{DynamicImage, RgbImage, GenericImageView, GenericImage, Rgb, FilterType};
use num_cpus;

const RESIZING_FACTOR: u32 = 1;
const ERROR_POWER: f32 = 2.0;
const SIZE_POWER: f32 = 0.0;
const AVG_INITAL_STEPS: f32 = 8.0;
const RELATIVE_MARGIN: f32 = 2.0;

struct Segment {
    pub image: DynamicImage,
    pub position: (i32, i32)
}

impl Segment {
    pub fn get_pixel(&self, x: i32, y: i32) -> Option<Rgb<u8>> {
        let position = ((x - self.position.0), (y - self.position.1));
        if position.0 >= 0 && position.1 >= 0 && position.0 < self.image.width() as i32 && position.1 < self.image.height() as i32 {
            match &self.image {
                DynamicImage::ImageRgb8(image) => {
                    Some(*image.get_pixel(position.0 as u32, position.1 as u32))
                },
                _ => None
            }
        } else {
            None
        }
    }
}

// Joins two segments and produces a new one
// O(n^2), single-threaded
fn join_segments(a: Segment, b: Segment) -> Segment {
    println!("Joining images");

    let start: (i32, i32) = (
        cmp::min(a.position.0, b.position.0),
        cmp::min(a.position.1, b.position.1)
    );
    let end: (i32, i32) = (
        cmp::max(a.position.0 + a.image.width() as i32, b.position.0 + b.image.width() as i32),
        cmp::max(a.position.1 + a.image.height() as i32, b.position.1 + b.image.height() as i32)
    );
    let size: (u32, u32) = (
        (end.0 - start.0) as u32,
        (end.1 - start.1) as u32
    );

    let mut image = RgbImage::new(size.0, size.1);
    
    for x in start.0..end.0 {
        for y in start.1..end.1 {
            let mut pixel = Rgb([0, 0, 0]);

            let pixel_a = a.get_pixel(x, y);
            let pixel_b = b.get_pixel(x, y);

            if let Some(pixel_a) = pixel_a {
                if let Some(pixel_b) = pixel_b {
 
                    let a = (
                        (a.position.0 + a.image.dimensions().0 as i32 / 2),
                        (a.position.1 + a.image.dimensions().1 as i32 / 2)
                    );
                    let b = (
                        (b.position.0 + b.image.dimensions().0 as i32 / 2),
                        (b.position.1 + b.image.dimensions().1 as i32 / 2)
                    );
                    let p = (
                        x as i32,
                        y as i32
                    );
                    let a_p = (p.0 - a.0, p.1 - a.1);
                    let a_b = (b.0 - a.0, b.1 - a.1);
                    let a_b_len = ((a_b.0 * a_b.0 + a_b.1 * a_b.1) as f32).sqrt();
                    let a_b_norm = (
                        a_b.0 as f32 / a_b_len,
                        a_b.1 as f32 / a_b_len
                    );
                    let dot = a_p.0 as f32 * a_b_norm.0 + a_p.1 as f32 * a_b_norm.1;
                    let c = (
                        a.0 + (a_b_norm.0 * dot) as i32,
                        a.1 + (a_b_norm.1 * dot) as i32
                    );
                    let a_b_dist = (b.0 - a.0).abs() + (b.1 - a.1).abs();
                    let a_c_dist = (c.0 - a.0).abs() + (c.1 - a.1).abs();
                    let b_c_dist = (c.0 - b.0).abs() + (c.1 - b.1).abs();
                    let a_c_b_dist = a_c_dist + b_c_dist;

                    if a_c_b_dist > a_b_dist {
                        if a_c_dist > b_c_dist {
                            pixel = pixel_b;
                        } else {
                            pixel = pixel_a;
                        }
                    } else {
                        assert_eq!(a_c_b_dist, a_b_dist);
                        pixel[0] = ((pixel_a[0] as i32 * b_c_dist + pixel_b[0] as i32 * a_c_dist) / a_c_b_dist) as u8;
                        pixel[1] = ((pixel_a[1] as i32 * b_c_dist + pixel_b[1] as i32 * a_c_dist) / a_c_b_dist) as u8;
                        pixel[2] = ((pixel_a[2] as i32 * b_c_dist + pixel_b[2] as i32 * a_c_dist) / a_c_b_dist) as u8;
                    }

                }
            }

            if pixel_a.is_none() {
                if let Some(pixel_b) = pixel_b {
                    pixel = pixel_b;
                }
            }

            if pixel_b.is_none() {
                if let Some(pixel_a) = pixel_a {
                    pixel = pixel_a;
                }
            }

            image.put_pixel((x - start.0) as u32, (y - start.1) as u32, pixel);
        }
    }

    image.save("result.jpg").unwrap();

    Segment {
        image: DynamicImage::ImageRgb8(image),
        position: (0, 0)
    }
}

// Adjusts the brightness of segment b to match segment a
// O(n^2), single-threaded
fn adjust_brightness(a: &Segment, b: &mut Segment) {
    println!("Adjusting brightness");

    let start: (i32, i32) = (
        cmp::min(a.position.0, b.position.0),
        cmp::min(a.position.1, b.position.1)
    );
    let end: (i32, i32) = (
        cmp::max(a.position.0 + a.image.dimensions().0 as i32, b.position.0 + b.image.dimensions().0 as i32),
        cmp::max(a.position.1 + a.image.dimensions().1 as i32, b.position.1 + b.image.dimensions().1 as i32)
    );

    let mut factor: (f32, f32, f32) = (0.0, 0.0, 0.0);
    let mut count: (u32, u32, u32) = (0, 0, 0);
    
    for x in start.0..end.0 {
        for y in start.1..end.1 {
            let pixel_a = a.get_pixel(x, y);
            let pixel_b = b.get_pixel(x, y);

            if let Some(pixel_a) = pixel_a {
                if let Some(pixel_b) = pixel_b {
                    if pixel_b[0] != 0 {
                        factor.0 += pixel_a[0] as f32 - pixel_b[0] as f32;
                        count.0 += pixel_b[0] as u32;
                    }
                    if pixel_b[1] != 0 {
                        factor.1 += pixel_a[1] as f32 - pixel_b[1] as f32;
                        count.1 += pixel_b[0] as u32;
                    }

                    if pixel_b[2] != 0 {
                        factor.2 += pixel_a[2] as f32 - pixel_b[2] as f32;
                        count.2 += pixel_b[0] as u32;
                    }
                }
            }
        }
    }

    factor.0 /= count.0 as f32;
    factor.1 /= count.1 as f32;
    factor.2 /= count.2 as f32;

    println!("Brightness factor is {:?}", factor);

    for x in 0..b.image.width() {
        for y in 0..b.image.height() {
            let mut pixel = b.image.get_pixel(x, y);
            pixel[0] = (pixel[0] as f32 * (1.0 + factor.0)) as u8;
            pixel[1] = (pixel[1] as f32 * (1.0 + factor.1)) as u8;
            pixel[2] = (pixel[2] as f32 * (1.0 + factor.2)) as u8;

            b.image.put_pixel(x, y, pixel);
        }
    }
    
}

// finds the position of image b in image a
// O(log(n^2) n^2), multi-threaded
fn find_position(a: &DynamicImage, b: &DynamicImage) -> (i32, i32) {
    println!("Finding position of image");

    let cpus = num_cpus::get();

    let mut best = ((256.0 * 3.0) as f32).powf(4.0);
    let mut position = (0, 0);

    let margin_x = (b.width() as f32 / RELATIVE_MARGIN) as u32;
    let margin_y = (b.height() as f32 / RELATIVE_MARGIN) as u32;

    let mut l_bound = - (margin_x as i32);
    let mut r_bound = a.width() as i32 + margin_x as i32 - b.width() as i32;
    let mut u_bound = - (margin_y as i32);
    let mut d_bound = a.height() as i32 + margin_y as i32 - b.height() as i32;

    let max_granularity = ((a.width() + a.height() + b.width() + b.height()) as f32 / 4.0 / AVG_INITAL_STEPS).log(2.0) as u32;
    //let max_granularity = 0;

    //println!("{} {} {} {}", l_bound, r_bound, u_bound, d_bound);

    println!("Max granularity is {}", max_granularity);

    for j in 0..=max_granularity {
        let mut threads = Vec::new();
        let step = (2 as i32).pow(max_granularity - j);

        for i in 0..cpus {
            let from = (r_bound - l_bound) / cpus as i32 * i as i32 + l_bound;
            let to = (r_bound - l_bound) / cpus as i32 * (i as i32 + 1) + l_bound;

            let a = a.clone();
            let b = b.clone();

            threads.push(thread::spawn(move || {
                let mut best = std::f32::MAX;
                let mut position = (0, 0);

                for ix in (from / step)..(to / step) {
                    for iy in (u_bound / step)..(d_bound / step) {
                        let x = ix * step;
                        let y = iy * step;
            
                        let l = cmp::max(0, x);
                        let r = cmp::min(a.width() as i32, x + b.width() as i32);
                        let u = cmp::max(0, y);
                        let d = cmp::min(a.height() as i32, y + b.height() as i32);
            
                        let size = (r - l, d - u);
                        let position_a = (l, u);
                        let position_b = (l - x, u - y);
                        
                        let mut diff = 0.0;
                        let mut n: u32 = 0;
            
                        for u in 0..size.0 {
                            for v in 0..size.1 {
                                let pixel_a = a.get_pixel((position_a.0 + u as i32) as u32, (position_a.1 + v as i32) as u32);
                                let pixel_b = b.get_pixel((position_b.0 + u as i32) as u32, (position_b.1 + v as i32) as u32);
                                let error = (
                                    ((pixel_a[0] as i32 - pixel_b[0] as i32).pow(2) + 
                                    (pixel_a[1] as i32 - pixel_b[1] as i32).pow(2) + 
                                    (pixel_a[2] as i32 - pixel_b[2] as i32).pow(2)) as f32
                                ).sqrt();
                                diff += error.powf(ERROR_POWER);
                                n += 1;
                            }
                        }
            
                        diff /= n as f32;
                        diff /= ((size.0 * size.1) as f32).powf(SIZE_POWER);
            
                        if diff < best {
                            best = diff;
                            position = (x, y);
                        }
                    }
                }

                (best, position)
            }));
        }

        for thread in threads {
            let (best_thread, position_thread) = thread.join().unwrap();

            if best_thread < best {
                best = best_thread;
                position = position_thread;
            }
        }

        println!("New best position with step size {} is {:?}", step, position);

        l_bound = cmp::max(- (margin_x as i32), position.0 - step as i32 - cpus as i32 / 2);
        r_bound = cmp::min(a.width() as i32 + margin_x as i32 - b.width() as i32, position.0 + step as i32 + cpus as i32 / 2);
        u_bound = cmp::max(- (margin_y as i32), position.1 - step as i32 - cpus as i32 / 2);
        d_bound = cmp::min(a.height() as i32 + margin_y as i32 - b.height() as i32, position.1 + step as i32 + cpus as i32 / 2);
    }

    position

    /*
    let mut l_bound = 0;
    let mut r_bound = a.width();
    let mut u_bound = 0;
    let mut d_bound = a.height();

    //let max_granularity = ((a.width() + a.height() + b.width() + b.height()) as f32 / 4.0 / AVG_INITAL_STEPS).log(2.0) as u32;
    let max_granularity = 0;

    println!("Max granularity is {}", max_granularity);

    for j in 0..=max_granularity {
        let mut threads = Vec::new();
        let granularity = (2 as u32).pow(max_granularity - j);

        for i in 0..cpus {
            let from = (r_bound - l_bound) / cpus as u32 * i as u32 + l_bound;
            let to = (r_bound - l_bound) / cpus as u32 * (i as u32 + 1) + l_bound;

            let a = a.clone();
            let b = b.clone();

            threads.push(thread::spawn(move || {
                let mut best = std::f32::MAX;
                let mut position = (0, 0);

                for ix in (from / granularity)..(to / granularity) {
                    let x = ix * granularity;
                    for iy in (u_bound / granularity)..(d_bound / granularity) {
                        let y = iy * granularity;
            
                        let l = cmp::max(b.width() / 2, x);
                        let r = cmp::min(b.width() / 2 + a.width(), x + b.width());
                        let u = cmp::max(b.height() / 2, y);
                        let d = cmp::min(b.height() / 2 + a.height(), y + b.height());
            
                        let size = (r - l, d - u);
                        let position_a = (l - b.width() / 2, u - b.height() / 2);
                        let position_b = (l - x, u - y);
                        
                        let mut diff = 0.0;
                        let mut n: u32 = 0;
            
                        for u in 0..size.0 {
                            for v in 0..size.1 {
                                let pixel_a = a.get_pixel(position_a.0 + u, position_a.1 + v);
                                let pixel_b = b.get_pixel(position_b.0 + u, position_b.1 + v);
                                let error = ((pixel_a[0] as i32 - pixel_b[0] as i32).abs() + (pixel_a[1] as i32 - pixel_b[1] as i32).abs() + (pixel_a[2] as i32 - pixel_b[2] as i32).abs()) as f32;
                                diff += error.powf(ERROR_POWER);
                                n += 1;
                            }
                        }
            
                        diff /= n as f32;
                        diff /= ((size.0 * size.1) as f32).powf(SIZE_POWER);
            
                        if diff < best {
                            best = diff;
                            position = (x as i32 - b.width() as i32 / 2, y as i32 - b.height() as i32 / 2);
                        }
                    }
                }

                (best, position)
            }));
        }

        for thread in threads {
            let (best_thread, position_thread) = thread.join().unwrap();

            if best_thread < best {
                best = best_thread;
                position = position_thread;
            }
        }

        println!("New best position with step size {} is {:?}", granularity, position);

        l_bound = cmp::max(0, (position.0 + b.width() as i32 / 2) - granularity as i32 - cpus as i32 / 2) as u32;
        r_bound = cmp::min(a.width() as i32, (position.0 + b.width() as i32 / 2) + granularity as i32 + cpus as i32 / 2) as u32;
        u_bound = cmp::max(0, (position.1 + b.height() as i32 / 2) - granularity as i32 - cpus as i32 / 2) as u32;
        d_bound = cmp::min(a.height() as i32, (position.1 + b.height() as i32 / 2) + granularity as i32 + cpus as i32 / 2) as u32;
    }
    
    position
    */
}

fn main() -> io::Result<()> {
    let mut result: Option<Segment> = None;

    println!("Number of virtual cores is {}", num_cpus::get());

    for entry in fs::read_dir("images")? {
        let path = entry?.path();
        println!("Loading image {}", path.to_str().unwrap());

        let image = image::open(path).unwrap();
        let image = image.resize(image.width() / RESIZING_FACTOR, image.height() / RESIZING_FACTOR, FilterType::Gaussian);

        let position = if let Some(r) = &result {
            find_position(&r.image, &image)
        } else {
            (0, 0)
        };
        println!("Position of image is {:?}", position);

        let mut segment = Segment {
            image,
            position
        };

        if let Some(r) = &result {
            adjust_brightness(&r, &mut segment);
        }

        result = if let Some(r) = result {
            Some(join_segments(r, segment))
        } else {
            Some(segment)
        };

    }

    println!("Saving result");

    result.unwrap().image.save("result.jpg")?;

    Ok(())
}
