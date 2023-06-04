/*
	Mandelbrot project 
	Hamish Hill <2002431@uad.ac.uk>
	Adapted from
		|Mandelbrot set example
		|Adam Sampson <a.sampson@abertay.ac.uk>
 */

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include <thread>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <condition_variable>  
#include <array>
#include <vector>

#define THREADS 12
#define THREADSM 100 //thread multiplier - determines total number of sections the mandlebrot set will be divided into

// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::complex;
using std::cout;
using std::endl;
using std::ofstream;
using std::thread;
using std::condition_variable;
using std::mutex;
using std::unique_lock;
// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_clock;

static std::atomic<int> aInt;
bool alive = true;
mutex m;

condition_variable task;
condition_variable complete;

// The size of the image to generate.
const int WIDTH = 1920;
const int HEIGHT = 1200;

// The number of times to iterate before we assume that a point isn't in the
// Mandelbrot set.
// (You may need to turn this up if you zoom further into the set.)
const int MAX_ITERATIONS = 500;

// The image data.
// Each pixel is represented as 0xRRGGBB.
uint32_t image[HEIGHT][WIDTH];

int waiting = 0;

// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char *filename)
{
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char *)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}

struct wThread {
	double left = -2.0;
	double right = 1.0;
	double top = 1.125;
	double bottom = -1.125;
	double length = HEIGHT / (THREADS*THREADSM);
	double startY = 0;
	double endY = 0;
};

void Calculate(struct wThread* args)
{
	for (int y = args->startY; y < args->endY; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			complex<double> c(args->left + (x * (args->right - args->left) / WIDTH),
				args->top + (y * (args->bottom - args->top) / HEIGHT));

			// Start off z at (0, 0).
			complex<double> z(0.0, 0.0);

			// Iterate z = z^2 + c until z moves more than 2 units
			// away from (0, 0), or we've iterated too many times.
			int iterations = 0;
			while (abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = (z * z) + c;

				++iterations;
			}

			if (iterations == MAX_ITERATIONS)
			{
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = 30; //0x000000; // black
			}
			else
			{
				// z escaped within less than MAX_ITERATIONS
				// iterations. This point isn't in the set.
				image[y][x] = 0xFFFF0000; // red

			}
		}
	}

}

int waitThread()
{
	unique_lock<mutex> lock(m);
	while (!waiting && alive) //no jobs but program still alive
	{
		task.wait(lock);	//thread waits for a job by default
	}
	if (!alive)
	{
		return -1;
	}

	else
	{
		return --waiting;
	}
}

void myThreadFunc1( std::vector<wThread>& args)
{
	int j;
	
	while (true)
	{

		j = waitThread();

		if (j == -1)
		{
			return;
		}

		//draw pixels
		Calculate(&args[j]);

		if (!--aInt)	//if unprocessed jobs reaches 0 
		{
			complete.notify_one();
		}
	}
}

void startWork(int &waiting, int jobs)
{
	unique_lock<mutex> lock(m);
	waiting = aInt = jobs;
	task.notify_all();			//start all threads
}

void waitForFinish()
{
	unique_lock<mutex> lock(m);
	while (aInt)		//wait for threads to finish work 
	{
		complete.wait(lock);
		
	}
	//cout << "threads complete" << endl;
}

void Finish()
{
	alive = false;
	task.notify_all(); //wake all threads for termination
}

int main(int argc, char *argv[])
{
	// Start timing
	the_clock::time_point start = the_clock::now();
	
	alive = true;

	std::array <thread, THREADS> myThread;
	std::vector <wThread> args;

	args.resize(THREADS*THREADSM);

	int jobs = THREADS * THREADSM; //total number of jobs

	aInt = 0;	//atomic int tracking number of incomplete jobs
	
	for (int i = 0; i < jobs; i++)	//organise list of jobs
	{
		double j = i;
		args.push_back(wThread());
		args[i].startY = j * args[i].length;
		args[i].endY = args[i].startY + args[i].length;
	}

	for(int i = 0; i < THREADS; i++)	//initiliase threads
	{
		myThread[i] = thread(myThreadFunc1, args);
	}

	//start work in all threads
	startWork(waiting, jobs);

	//main waits until no work is left for threads
	waitForFinish();

	//when no work is left threads are woken for termination
	Finish();

	for (int i = 0; i <THREADS; i++)
	{
		
		myThread[i].join();
	}

	// Stop timing
	the_clock::time_point end = the_clock::now();

	// Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Computing the Mandelbrot set took " << time_taken << " ms." << endl;

	write_tga("output.tga");
	return 0;
}
