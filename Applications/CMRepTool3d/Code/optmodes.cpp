#include "glui.h"

// Optimization driver
SplineOptDriver *optimizationDriver = NULL;
SplineDistanceMatcher *distanceMatcher = NULL;

void OptimizationMH::onDraw() {
	// Call the print command
	// executeCommand("match");

	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	glColor3d(0.75,0.75,0.75);
	glRasterPos2d(GLDisplayDriver::width-500,25);	
	fntCourier18->printf("Patch %2d,%2d  DIST = %7lg  EVALS = %7lg",u,v,value,cpStartOptCost+optCost);
	glPopAttrib ();
}

OptimizationMH::OptimizationMH() {
	window = 0;
	msPerWindow = 1000;
}

void OptimizationMH::makeControlSequence() {
	// Create a sequence of control point indexes
	cpList.clear();
	for(int v=0;v<=spline->dim(1)-window+1;v++) {
		int base = v << 16;
		for(int u=0;u<=spline->dim(0)-window+1;u++) {
			cpList.push_back(base + u);
		}
	}
	
	// Shuffle the control points
	random_shuffle(cpList.begin(),cpList.end());
}

void OptimizationMH::start() {
	// Clear all parameters for idler
	cpTimeLeft = 0;

	// Create needed objects
	od = NULL;
	match = imaging.getImageMatch();

	// Bail if there is no match
	if(!match) {
		throw "Image is not loaded or is of incorrect type";
	}
	match->setSpline(spline,splineData);

	// Save the total value
	value = match->getMatch();
	cpStartOptCost = 0;
	optCost = 0;
	cpTimeLeft = 0;

	// Get the window size
	window = settings->getIntValue("optimizer.windowSize",0);

	// Add an idle display listener
	GLDisplayDriver::addListener(this,GLDisplayDriver::IDLE);	
	GLDisplayDriver::addRenderer(this,GLDisplayDriver::EYE);
}

void OptimizationMH::stop() {
	// Delete the optimization driver and the distance map
	if(od) delete od;
	od = NULL;

	// Remove the idle display listener
	GLDisplayDriver::removeRenderer(this,GLDisplayDriver::EYE);
	GLDisplayDriver::removeListener(this,GLDisplayDriver::IDLE);	
}

bool OptimizationMH::handleIdle() {
	if(cpTimeLeft <= 0) {
		
		// Delete old optimization driver
		if(od) {
			delete od;
		}
		
		// Update total cost
		cpStartOptCost += optCost;

		// If we are doing windowed optimization, pick a window
		if(window > 0) {
			// If there is no more time for current control point, select next control point
			if(cpList.size() == 0) {
				makeControlSequence();
				msPerWindow += 250;
			}

			// Pop off from the list
			int uv = cpList.back();
			cpList.pop_back();
			u = uv & 0xffff;
			v = uv >> 16;

			// Create an optimization driver
			od = new SplineOptDriver(spline,splineData,match,u,v,window,window,&settings->getSubFolder("optimizer"));
		}
		else {
			u = 0;
			v = 0;
			od = new SplineOptDriver(spline,splineData,match,0,0,spline->dim(0)+1,spline->dim(1)+1,
				&settings->getSubFolder("optimizer"));
		}

		// Assign a new wait time
		cpTimeLeft = msPerWindow;
	}

	// Perform an optimization
	if(od->optimize(250)) {
		cpTimeLeft = 0;
	} else {
		cpTimeLeft -= 250;
	}

	// Apply the best solution
	od->applyBestSolution();

	// Store the value
	value = od->getBestValue();
	optCost = od->getEvalCost();

	// Do not block other idlers
	return false;
}


void RigidMatchMode::onDraw() {
	if(match && es) {		
		float lkhd = match->getMatch();
		float post = es->getBestEverValue();
		double cost = srp->getEvaluationCost();

		glPushAttrib(GL_LIGHTING_BIT);
		glDisable(GL_LIGHTING);
		glColor3d(0.75,0.75,0.75);
		glRasterPos2d(GLDisplayDriver::width-500,25);	
		fntCourier18->printf("LKHD = %7g  POST = %7g  EVALS = %7lg",lkhd,post,cost);
		glPopAttrib ();
	}
}

void RigidMatchMode::start() {
	match = imaging.getImageMatch();

	// Bail if there is no match
	if(!match) {
		throw "Image is not loaded or is of incorrect type";
	}
	match->setSpline(spline,splineData);

	srp = new SplineRigidProblem(spline,splineData,match);

	Vector mean(7),sd(7);
	srp->makeVector(mean);
	srp->makeSDVector(sd);
	sd = sd * 0.1;

	ss = new GaussianSS(mean,sd);
	es = new EvolutionaryStrategy(*srp,*ss,
		settings->getIntValue("optimizer.es.mu",2),
		settings->getIntValue("optimizer.es.lambda",4),SELECTION_MuPlusLambda);
	es->setNth(0,mean);

	// Set upper/lower bounds
	Vector lb(7,-3.15,-3.15,-3.15,-1.0,-1.0,-1.0,-1.0);
	Vector ub(7,3.15,3.15,3.15,1.0,1.0,1.0,1.0);
	es->setXBounds(lb,ub);

	// Add an idle display listener
	GLDisplayDriver::addRenderer(this,GLDisplayDriver::EYE);
	GLDisplayDriver::addListener(this,GLDisplayDriver::IDLE);	
}

void RigidMatchMode::stop() {
	delete srp;
	delete ss;
	delete es;

	// Remove the idle display listener
	GLDisplayDriver::removeRenderer(this,GLDisplayDriver::EYE);
	GLDisplayDriver::removeListener(this,GLDisplayDriver::IDLE);	
}

bool RigidMatchMode::handleIdle() {
	es->performIteration();
	srp->applyVector(es->getBestEverX());
	return false;
}
/*
void CurveAxisMatchMode::start() {

}

void CurveAxisMatchMode::stop() {

}

bool CurveAxisMatchMode::handleIdle() {
	return false;
}

void CurveAxisMatchMode::onDraw() {

}*/
