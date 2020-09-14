use rand::Rng;
use std::f64::consts::PI;
use num::complex::Complex;
use plotlib::page::Page;
use plotlib::repr::Plot;
use plotlib::view::ContinuousView;
use plotlib::style::{PointMarker, PointStyle};

static STEPS: usize = 100;
static COLS: usize = 700;
static PERIODS: usize = 1024;

fn linspace(start: f64, stop: f64, step: usize) -> Vec<f64> {
	let mut res_vec: Vec<f64> = vec![start];
	let size = (stop-start)/(step as f64);

	for i in 1..step+1 {
		res_vec.push(start + (size*(i as f64)));
	}
	res_vec
}

#[derive(Debug)]
struct ECA<'a> {
	n: &'a usize,
	t: &'a usize,
	rule: &'a Vec<usize>
}

#[derive(Debug)]
struct Spectrum<'b> {
	n: &'b usize,
	t: &'b usize,
	x: &'b Vec<Vec<usize>>,
	freq: &'b Vec<f64>
}

fn eca_rule_result(x: (usize, usize, usize), rule: Vec<usize>) -> usize {
	let rule_loc = match x {
		(0, 0, 0) => 7,
		(0, 0, 1) => 6,
		(0, 1, 0) => 5,
		(0, 1, 1) => 4,
		(1, 0, 0) => 3,
		(1, 0, 1) => 2,
		(1, 1, 0) => 1,
		(1, 1, 1) => 0,
		_ => panic!("Not allowed"),
	};
	if rule[rule_loc] == 1 {
		1
	} else {
		0
	}
}

fn construct(n: &usize, t: &usize, rule: &Vec<usize>) -> Vec<Vec<usize>> {
	let mut x = vec![vec![0; *n]; *t];
	let mut rng = rand::thread_rng();

	for i in 0..*n {
		x[0][i] = rng.gen_range(0,2);
	}
	for m in 1..*t {
		for i in 1..*n-1 {
			x[m][i] = eca_rule_result((x[m-1][i-1],x[m-1][i],x[m-1][i+1]), rule.to_vec());
		}
	}
	x
}

fn dft(n: usize, t: &usize, x: &Vec<Vec<usize>>, f: f64) -> Complex<f64> {
	let mut x_hat = Complex::new(0.0,0.0);
	for j in 0..*t {
		x_hat.re += (x[j][n] as f64)*(2.0*PI*(*t as f64)*f/(*t as f64)).cos()/(*t as f64);
		x_hat.im += (x[j][n] as f64)*(2.0*PI*(*t as f64)*f/(*t as f64)).sin()/(*t as f64);
	}
	x_hat
}

fn density (dft_vec: Vec<Vec<Complex<f64>>>) -> Vec<f64> {
	let mut s = vec![0.0; STEPS];
	for f in 0..STEPS {
		for i in 0..COLS{
			s[f] += (dft_vec[f][i].re*dft_vec[f][i].re)+(dft_vec[f][i].im*dft_vec[f][i].im);
		}
		s[f] = s[f]/(COLS as f64);
	}
	s
}

fn main() {
	let rule = vec![0, 1, 1, 0, 1, 1, 1, 0];
	let n = COLS;
	let t = PERIODS;

	let x = ECA{ n: &n, t: &t, rule: &rule };
	let y = construct(x.n,x.t,x.rule);

	let freq = linspace(0.0,10.0,STEPS);
	
	let spectral = Spectrum { n: &n, t: &t, x: &y, freq: &freq} ;

	let k = &n;
	let mut dft_vec = vec![vec![Complex::new(0.0,0.0); *k]; freq.len()];

	for j in 0..freq.len() {
		for i in 0..n {
			dft_vec[j][i] = dft(i, spectral.t, spectral.x, spectral.freq[j]);
		}
	}

	let s = density(dft_vec);
	let mut res_vec: Vec<(f64,f64)> = vec![((freq[0],s[0]))];

	for i in 1..STEPS {
		res_vec.push((freq[i],s[i]));
	}

	let data1 = res_vec;

    // We create our scatter plot from the data
    let s1: Plot = Plot::new(data1).point_style(
        PointStyle::new()
            .marker(PointMarker::Square) // setting the marker to be a square
            .colour("#DD3355"),
    ); // and a custom colour

    // The 'view' describes what set of data is drawn
    let v = ContinuousView::new()
        .add(s1)
        .x_range(0., 10.)
        .y_range(0.330, 0.334)
        .x_label("Frequency")
        .y_label("Power Spectrum");

    // A page with a single view is then saved to an SVG file
    Page::single(&v).save("scatter.svg").unwrap();
}
