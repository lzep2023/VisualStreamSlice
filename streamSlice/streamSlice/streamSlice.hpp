#ifndef STREAM_SLICE_HPP
#define STREAM_SLICE_HPP
#include"def.h"
#define eps 2.2204e-16

std::vector<cv::Point2f> traceStream(cv::Mat1f xvec, cv::Mat1f yvec, cv::Mat1f ugrid, cv::Mat1f vgrid, int xdim, int ydim, float sx, float sy, float step, int maxVert)
{
	int numVerts = 0;
	float x = sx;
	float y = sy;
	int ix, iy;
	float x0, x1, y0, y1, xi, yi, dx, dy;
	float xfrac, yfrac, ui, vi;
	float a, b, c, d, imax;
	std::vector<cv::Point2f> verts(0);
	while (1)
	{
		if (x < 0 || x>xdim - 1 || y<0 || y>ydim - 1 || numVerts >= maxVert)
			break;

		if (std::isnan(x))
		{
			ix = 0;
		}
		else
		{
			ix = x;
		}

		if (std::isnan(x))
		{
			iy = 0;
		}
		else
		{
			iy = y;
		}
		/*if (ix < 0) 
		{
			ix = 0; x = 1;
		}
		if (iy < 0)
		{
			iy = 0;
			y = 1;
		}*/
		if (ix >= xdim - 1) { x = xdim - 1; ix = x-1; }
		if (iy >= ydim - 1) {
			y = ydim - 1;
			iy=y-1;
		}
		xfrac = x - ix;
		yfrac = y - iy;

		a = (1 - xfrac)*(1 - yfrac);
		b = (xfrac)*(1 - yfrac);
		c = (1 - xfrac)*(yfrac);
		d = (xfrac)*(yfrac);

		/*x0 = xvec(iy, ix);
		x1 = xvec(iy, ix + 1);
		y0 = yvec(iy, ix);
		y1 = yvec(iy + 1, ix);*/
		x0 = ix;
		x1 = ix + 1;
		y0 = iy;
		y1 = iy+1;
		xi = x0*(1 - xfrac) + x1*xfrac;
		yi = y0*(1 - yfrac) + y1*yfrac;

		verts.push_back(cv::Point2f(xi,yi));  //x+1, y+1

		int sz = verts.size();
		numVerts = sz;
		//std::cout << "numVerts:" << numVerts << std::endl;
		if (sz > 2)
		{
			float m1 = verts[sz-1].x;
			float n1 = verts[sz-1].y;

			float m2 = verts[sz-3].x;
			float n2 = verts[sz-3].y;

			float diffm12 = fabs(m1 - m2);
			float diffn12 = fabs(n1 - n2);

			if (diffm12==0&&diffn12==0)
				break;
		}

		/*if (sz >= 500)
		{
			break;
		}*/

	 
		ui = ugrid(iy, ix)*a + ugrid(iy, ix + 1)*b + ugrid(iy + 1,ix)*c + ugrid(iy + 1, ix + 1)*d;
		vi = vgrid(iy, ix)*a + vgrid(iy, ix + 1)*b + vgrid(iy + 1,ix)*c + vgrid(iy + 1, ix + 1)*d;

		dx = x1 - x0;
		dy = y1 - y0;
		if (dx) ui /= dx;
		if (dy) vi /= dy;

		if (fabs(ui)>fabs(vi)) imax = fabs(ui); else imax = fabs(vi);
		if (imax == 0) break;

		imax = step / imax;

		ui *= imax;
		vi *= imax;

		x += ui;
		y += vi;
	}

	return verts;
}
std::vector<cv::Point2f> stream2(cv::Mat1f x, cv::Mat1f y, cv::Mat1f u, cv::Mat1f v, float sx, float sy, float step=0.1f, int maxVert=10000)
{
	return traceStream(x, y, u, v, u.cols, u.rows, sx, sy, step, maxVert);
}
void niceStreams(cv::Mat1f x, cv::Mat1f y, cv::Mat1f u, cv::Mat1f v, f32 density, bool arrows, std::vector< std::vector<cv::Point2f> >& linePt, std::vector< std::vector<cv::Point2f> >& arrowPt, int& gszx, int& gszy)
{
	int sv = std::min(x.rows, x.cols);
	float stepsize = std::min(0.1, ((sv - 1)*1.0) / 100);
	int maxVert = 10000;

	int num = 20;

	int nd = ceil(num*density);

	int nrstart = nd;
	int ncstart = nd;
	int nrend = ceil(num*density * 4);
	int ncend = nrend;
	int nrarrow = nd;
	int ncarrow = nd;
	double xmin0, xmax0, ymin0, ymax0;
	float xmin, xmax, ymin, ymax, xrange, yrange;
	cv::minMaxIdx(x, &xmin0, &xmax0);
	cv::minMaxIdx(y, &ymin0, &ymax0);
	
	xmin = xmin0;
	ymin = ymin0;
	xmax = xmax0;
	ymax = ymax0;

	xrange = xmax - xmin;
	yrange = ymax - ymin;

	

	float incsx = xrange / ncstart;  //增长步长
	float incsy = yrange / nrstart;
	float ixrgs = ncstart / xrange*(1 - eps);
	float iyrgs = nrstart / yrange*(1 - eps);
	float ixrge = ncend / xrange*(1 - eps);
	float iyrge = nrend / yrange*(1 - eps);


	uMat startgrid = uMat::zeros(nrstart, ncstart);
	uMat endgrid = uMat::zeros(nrend, ncend);

	std::cout << "nrstart:" << nrstart << ",ncstart:" << ncstart << std::endl;
	std::cout << "nrend:" << nrend << ",ncend:" << ncend << std::endl;
	 

	uMat arrowgrid;
	float arrowscale;
	float ixrca;
	float iyrca;
	if (arrows)
	{
		arrowgrid = uMat::ones(nrarrow, ncarrow);
		for (int i = 1; i < nrarrow;i+=3)
		for (int j = 1; j < ncarrow; j += 3)
		{
			arrowgrid(i, j) = 0;
		}
		arrowscale = 0.01f*(xrange + yrange) / density;
		ixrca = (ncarrow*1.f) / xrange*(1 - eps);
		iyrca = (nrarrow*1.f) / yrange*(1 - eps);
	}

	int row = startgrid.rows;
	int col = startgrid.cols;
	for (int r = 0; r < row;r++)
	for (int c = 0; c < col; c++)
	{
		std::cout << "r:" << r << ",c:" << c << endl;
		if (startgrid(r, c) == 0)
		{
			startgrid(r, c) = 1;

			//float xstart = xmin + (c - 0.5)*incsx;
			//float ystart = ymin + (r - 0.5)*incsy;

			float xstart = xmin + (c - 0)*incsx+0.5;
			float ystart = ymin + (r - 0)*incsy+0.5;

			std::cout << "f1" << std ::endl;
			std::vector<cv::Point2f> vertsf = stream2(x, y, u, v, xstart, ystart,stepsize,maxVert);
			std::cout << "f2" << std::endl;
			std::vector<cv::Point2f> vertsb = stream2(x, y, -u, -v, xstart, ystart, stepsize, maxVert);
			std::cout << "f3" << std::endl;
			
			std::vector<cv::Point2f> vertsf0(0);
			for (int f = 0; f < vertsf.size(); f++)
			{
				if (std::isnan(vertsf[f].x) || std::isnan(vertsf[f].y))
					continue;
				else
					vertsf0.emplace_back(vertsf[f]);
			}
			//std::cout << "vertsf0 sz:" << vertsf0.size() << std::endl;
			std::vector<cv::Point2f> vertsb0(0);
			for (int b = 0; b < vertsb.size(); b++)
			{
				if (std::isnan(vertsb[b].x) || std::isnan(vertsb[b].y))
					continue;
				else
					vertsb0.emplace_back(vertsb[b]);
			}
			//std::cout << "vertsb0 sz:" << vertsf0.size() << std::endl;
			//forward
			std::vector<cv::Point2f> flpt(0);
			int tcc, trr,cc,rr;
			if (vertsf0.size()>0)
			{
				tcc = (int)(floorf((vertsf0[0].x - xmin)*ixrge));
				trr = (int)(floorf((vertsf0[0].y - ymin)*iyrge));
			}
			int p;
			for (p = 0; p < vertsf0.size(); p++)
			{
				float xc = vertsf0[p].x;
				float yc = vertsf0[p].y;

				cc = floorf((xc - xmin)*ixrgs);
				rr = floorf((yc - ymin)*iyrgs);

				cc = (std::min)(cc, ncstart - 1);
				rr = (std::min)(rr, nrstart - 1);
				startgrid(rr, cc) = 1; //设置grid为1

				cc = floor((xc - xmin)*ixrge);
				rr = floor((yc - ymin)*iyrge);
				cc = (std::min)(cc, ncend - 1);
				rr = (std::min)(rr, nrend - 1);
				if (endgrid(rr, cc) == 1)
				{
					if (cc!=tcc||rr!=trr)
					{
						break;
					}
				}
				else
				{
					tcc = cc;
					trr = rr;
					endgrid(rr, cc) = 1; //设置endgrid
				}

				if (arrows)
				{
					cc = floorf((xc - xmin)*ixrca);
					rr = floorf((yc - ymin)*iyrca);
					cc = (std::min)(cc, ncarrow - 1);
					rr = (std::min)(rr, nrarrow - 1);

					float dx,dy,dnorm;
					if (p > 0 && arrowgrid(rr, cc) == 0)
					{
						arrowgrid(rr, cc) = 1;
						dx = vertsf0[p].x-vertsf0[p - 1].x;
						dy = vertsf0[p].y-vertsf0[p - 1].y;
						dnorm = sqrtf(dx*dx + dy*dy)+eps;
						dx /= dnorm;
						dy /= dnorm;

						float arrowlen = 1;
						float arrowwidth = 0.6f;
						float p2x = xc - arrowlen*dx*arrowscale;
						float p2y = yc - arrowlen*dy*arrowscale;

						float p2xprev = p2x - arrowwidth*(-dy)*arrowscale;
						float p2yprev = p2y - arrowwidth*(dx)*arrowscale;
						float p2xbef = p2x + arrowwidth*(-dy)*arrowscale;
						float p2ybef = p2y + arrowwidth*(dx)*arrowscale;

						//将箭头点存档
						std::vector<cv::Point2f> arrpt(0);
						arrpt.emplace_back(cv::Point2f(p2xprev, p2yprev));
						arrpt.emplace_back(cv::Point2f(xc, yc));
						arrpt.emplace_back(cv::Point2f(p2xbef, p2ybef));
						arrowPt.push_back(arrpt);
					}
				}
			}

			std::cout << "copy" << std::endl;
			if (p>1)
			flpt.assign(vertsf0.begin(), vertsf0.begin() + p - 2);
			//backward
			std::vector<cv::Point2f> blpt(0);
			if (vertsb0.size()>0)
			{
				tcc = (int)(floorf((vertsb0[0].x - xmin)*ixrge));
				trr = (int)(floorf((vertsb0[0].y - ymin)*iyrge));
			}
			std::cout << "copy1" << std::endl;
			for (p = 0; p < vertsb0.size(); p++)
			{
				float xc = vertsb0[p].x;
				float yc = vertsb0[p].y;

				cc = floorf((xc - xmin)*ixrgs);
				rr = floorf((yc - ymin)*iyrgs);
				cc = (std::min)(cc, ncstart - 1);
				rr = (std::min)(rr, nrstart - 1);

				startgrid(rr, cc) = 1; //设置grid为1

				cc = floor((xc - xmin)*ixrge);
				rr = floor((yc - ymin)*iyrge);
				cc = (std::min)(cc, ncend - 1);
				rr = (std::min)(rr, nrend - 1);
				if (endgrid(rr, cc) == 1)
				{
					if (cc != tcc || rr != trr)
					{
						break;
					}
				}
				else
				{
					tcc = cc;
					trr = rr;
					endgrid(rr, cc) = 1; //设置endgrid
				}

				if (arrows)
				{
					cc = floorf((xc - xmin)*ixrca);
					rr = floorf((yc - ymin)*iyrca);
					cc = (std::min)(cc, ncarrow - 1);
					rr = (std::min)(rr, nrarrow - 1);

					float dx, dy, dnorm;
					if (p > 0 && arrowgrid(rr, cc) == 0)
					{
						arrowgrid(rr, cc) = 1;
						dx = vertsb0[p].x - vertsb0[p - 1].x;
						dy = vertsb0[p].y - vertsb0[p - 1].y;
						dx = -dx;
						dy = -dy;
						dnorm = sqrtf(dx*dx + dy*dy)+eps;
						dx /= dnorm;
						dy /= dnorm;

						float arrowlen = 1;
						float arrowwidth = 0.6f;
						float p2x = xc - arrowlen*dx*arrowscale;
						float p2y = yc - arrowlen*dy*arrowscale;

						float p2xprev = p2x - arrowwidth*(-dy)*arrowscale;
						float p2yprev = p2y - arrowwidth*(dx)*arrowscale;
						float p2xbef = p2x + arrowwidth*(-dy)*arrowscale;
						float p2ybef = p2y + arrowwidth*(dx)*arrowscale;

						//将箭头点存档
						std::vector<cv::Point2f> arrpt(0);
						arrpt.emplace_back(cv::Point2f(p2xprev, p2yprev));
						arrpt.emplace_back(cv::Point2f(xc, yc));
						arrpt.emplace_back(cv::Point2f(p2xbef, p2ybef));
						arrowPt.push_back(arrpt);
					}
				}
			}
			std::cout << "p:" << p << std::endl;
			for (int m = 0; m < p; m++)
			{
				blpt.emplace_back(vertsb0[p - m - 1]);
			}


		
			if (blpt.size()>0)
			linePt.push_back(blpt);
			if (flpt.size()>0)
			linePt.push_back(flpt);
		}
	}
}
void streamSlice2D(cv::Mat1f x, cv::Mat1f y, cv::Mat1f u, cv::Mat1f v, f32 density, bool arrows, cv::Mat3b& slImg)
{
	f32 _density = sqrtf(density); //线条之间的密度
	std::vector< std::vector<cv::Point2f> > linePt;
	std::vector< std::vector<cv::Point2f> > arrowPt;

	int gszx, gszy;
	std::cout << "start nice stream" << std::endl;
	niceStreams(x, y, u, v, _density, arrows, linePt, arrowPt, gszx, gszy);
	std::cout << "finish nice stream" << std::endl;
	float scale = 640.f / u.cols;

	slImg.create(u.rows*scale, 640);
	slImg.setTo(cv::Scalar(255, 255, 255));
	std::cout << "linePt.size():" << linePt.size() << std::endl;
	for (int lp = 0; lp < linePt.size(); lp++)
	{
		for (int q = 1; q < linePt[lp].size(); q++)
		{
			cv::line(slImg, cv::Point2f(linePt[lp][q - 1].x*scale, linePt[lp][q - 1].y*scale), cv::Point2f(linePt[lp][q].x*scale, linePt[lp][q].y*scale), cv::Scalar(0, 0, 0), 1);
		}
	}
	std::cout << "arrowPt.size():" << arrowPt.size() << std::endl;
	for (int ap = 0; ap < arrowPt.size(); ap++)
	{
		//std::cout << "arrowPt[ap].size():" << arrowPt[ap].size() << std::endl;
		for (int q = 1; q < arrowPt[ap].size(); q++)
		{
			cv::line(slImg, cv::Point2f(arrowPt[ap][q - 1].x*scale, arrowPt[ap][q - 1].y*scale), cv::Point2f(arrowPt[ap][q].x*scale, arrowPt[ap][q].y*scale), cv::Scalar(0, 0, 255), 1);
		}
	}

}
#endif