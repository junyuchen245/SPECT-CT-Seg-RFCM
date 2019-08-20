#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "miputil.h"

typedef short int LABELTYPE;
#define INDEX3D(x,y,z, xdim, ydim, zdim) ((xdim)*((z)*(ydim) + (y)) + (x))

//////////////////////////////////////////////////////////////////////////// k means ////////////////////////////////////////////////////////////////////////////
float max(float* voxvals, int nvox)
{
	float max_val = voxvals[0];
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]>max_val)
		{
			max_val = voxvals[i];
		}
	}
	return max_val;
}

float min(float* voxvals, int nvox)
{
	float min_val = voxvals[0];
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]<min_val)
		{
			min_val = voxvals[i];
		}
	}
	return min_val;
}

int argmin(float* voxvals, int nvox)
{
	float min_val = voxvals[0];
	int min_i = 0;
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]<min_val)
		{
			min_val = voxvals[i];
			min_i = i;
		}
	}
	return min_i;
}

int argmin_double(double* voxvals, int nvox)
{
	double min_val = voxvals[0];
	int min_i = 0;
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]<min_val)
		{
			min_val = voxvals[i];
			min_i = i;
		}
	}
	return min_i;
}

int argmax(float* voxvals, int nvox)
{
	float max_val = voxvals[0];
	int max_i = 0;
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]>max_val)
		{
			max_val = voxvals[i];
			max_i = i;
		}
	}
	return max_i;
}

int argmax_double(double* voxvals, int nvox)
{
	double max_val = voxvals[0];
	int max_i = 0;
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]>max_val)
		{
			max_val = voxvals[i];
			max_i = i;
		}
	}
	return max_i;
}

double mean_class(float* voxvals, int nvox, LABELTYPE* label_map, LABELTYPE label)
{
	double sum = 0.0;
	double count = 0;
	for(int i=0; i<nvox; i++)
	{
		if (label_map[i] == label)
		{
			sum += voxvals[i];
			count += 1.0;
		}
	}
	if(count==0)
	return 0;
	return sum/count;
}

double var_class(float* voxvals, int nvox, LABELTYPE* label_map, LABELTYPE label)
{
	double sum = 0.0;
	double count = 0;
	double mean = 0.0;
	double var  = 0.0;
	double nval = 0.0;

	for(int i=0; i<nvox; i++)
	{
		if (label_map[i] == label)
		{
			sum += voxvals[i];
			count += 1.0;
		}
	}
	mean = sum/count;

	sum = 0.0;
	for(int i=0; i<nvox; ++i){
		if (label_map[i] == label){
			sum += pow(mean-voxvals[i],2.0);
		}
	}

	if(count == 0)
	{
		nval = 0.0001;
	}
	else
	{
		nval = (double)count;
	}
	var=sum/nval; //Use the population std dev formula
	//fprintf(stderr, "var is:%f\n", var);
	return var;	
}

void k_means(float* voxvals, LABELTYPE* labels, int xdim, int ydim, int zdim, int num_clus, float delta, int simpleKmeans)
{
	int nvox = xdim*ydim*zdim;
	float part_lb;
	float part_ub;
	part_lb = min(voxvals, nvox); // Max value within the image
	part_ub = max(voxvals, nvox); // Min value within the image
	float ini_means[num_clus];
	for(int i=0; i<num_clus; i++)
	{
		ini_means[i] = ((part_ub - part_lb)/num_clus *((i+1)-0.5) + part_lb);
	}
	
	float means[num_clus];
	//means = &ini_means;
	for(int i=0; i<num_clus; i++)
	{
		means[i] = ini_means[i];
	}
	//fprintf(stderr, "mean1=%f, mean2=%f\n", means[0], means[1]);
	LABELTYPE* label_map;
	LABELTYPE* new_map;
	label_map = (LABELTYPE *)malloc(xdim*ydim*zdim*sizeof(LABELTYPE));
	new_map   = (LABELTYPE *)malloc(xdim*ydim*zdim*sizeof(LABELTYPE));
	int first_enter_flag = 1;
	int max_iter = 40;
	float val = 0.0;
	float pixel_pot[num_clus];
	float penalty = 0.0;
	LABELTYPE label;

	for(int iter=0;iter<max_iter;iter++)    //iteration
	{
		for(int ix=0;ix<xdim;ix++)          //x dim
		{
			for(int iy=0;iy<ydim;iy++)      //y dim
			{
				for(int iz=0;iz<zdim;iz++)  //z dim
				{
					val=voxvals[INDEX3D(ix, iy, iz, xdim, ydim, zdim)];
					for(LABELTYPE k=0;k<num_clus;k++)
					{
						penalty = 0.0;
						if(first_enter_flag == 0 && simpleKmeans == 0)
						{
							for(int ia=-1; ia<2; ia++)
							{
								for(int ib=-1; ib<2; ib++)
								{
									for(int ic=-1; ic<2; ic++)
									{
										if(ix + ia < 0 || iy + ib < 0 || iz + ic < 0 || ix + ia >= xdim || iy + ib >= ydim || iz + ic >= zdim) continue;
										label = label_map[INDEX3D(ix+ia, iy+ib, iz+ic, xdim, ydim, zdim)];
										if(label != k)
										{
											penalty = penalty + 1.0;
										}
									}//ic
								}//ib
							}//ia
						}//if first enter
						float smooth_pot = delta * penalty;
						float ini_pot = (val-means[k])*(val-means[k]);
						pixel_pot[k] = ini_pot + smooth_pot;
					}//class
					new_map[INDEX3D(ix, iy, iz, xdim, ydim, zdim)] = argmin(pixel_pot, num_clus);
				}//z
			}//y
		}//x
		for(int i=0; i < xdim*ydim*zdim; i++) label_map[i] = new_map[i];

		// recompute mean
		for(LABELTYPE k=0; k<num_clus; k++)
		{
			means[k] = mean_class(voxvals, nvox, label_map, k);
		}
		//fprintf(stderr, "k means: mean1=%f, mean2=%f, mean3=%f\n", means[0], means[1], means[2], means[3]);
		first_enter_flag = 0;
	}//iter
	for(LABELTYPE k=0; k<num_clus; k++)
	{
		int max_i;
		max_i = argmax(means, num_clus);
		means[max_i] = -10000.0;
		for(int i=0; i<nvox; i++)
		{
			if(label_map[i]==max_i)
			{
				labels[i] = num_clus - k - 1;
			}
		}
	}
	
	free((void *)new_map);
	free((void *)label_map);
	return;
}
void dot_multiply(double* results, double* vec_1, double* vec_2, int dim)
{
	for(int i=0; i<dim; i++)
	{
		results[i] = vec_1[i] * vec_2[i];
	}
}
void dot_divide(double* results, double* vec_1, double* vec_2, int dim)
{
	for(int i=0; i<dim; i++)
	{
		results[i] = vec_1[i] / vec_2[i];
	}
}
void dot_minus(double* results, double* vec_1, double* vec_2, int dim)
{
	for(int i=0; i<dim; i++)
	{
		results[i] = vec_1[i] - vec_2[i];
	}
}

int FCM(float* voxvals, LABELTYPE* labels, int xdim, int ydim, int zdim, int num_clus, double beta, double q)
{
	LABELTYPE* pixel_reg;
	int nvox     = xdim*ydim*zdim; // number of voxels
	int num_neib = 26; // number of neighborhood pixels
	double* centroids;
	double* u_jk_pre;
	double* u_jk;
	double* u_lm;
	double* u_lm_tmp;
	double power = -1/(q-1);
	pixel_reg = (LABELTYPE *)malloc(xdim*ydim*zdim*sizeof(LABELTYPE)); // define label map
	centroids = (double *)dvector(num_clus*sizeof(double)); // class means
	u_jk_pre  = (double *)malloc(xdim*ydim*zdim*num_clus*sizeof(double));
	u_jk      = (double *)malloc(xdim*ydim*zdim*num_clus*sizeof(double));
	u_lm      = (double *)malloc(num_clus*num_neib*sizeof(double));
	u_lm_tmp  = (double *)malloc(num_clus*num_neib*sizeof(double));

	for(int i=0; i<xdim*ydim*zdim*num_clus; i++)
	{
		u_jk[i]     = 0;
		u_jk_pre[i] = 0;
	}
	//fprintf(stderr, "class num is:%d\n", num_clus);
	k_means(voxvals, pixel_reg, xdim, ydim, zdim, num_clus, 0, 1); // k means clustering as initialization
	//for(int i=0; i < xdim*ydim*zdim; ++i)
	//	labels[i] = label_map[i];
	//return 0;

	for(int k=0; k<num_clus; k++)
	{
		centroids[k]     = mean_class(voxvals, nvox, pixel_reg, k);
	}//class k
	//fprintf(stderr, "initial: mean1=%f, mean2=%f, mean3=%f\n", centroids[0], centroids[1], centroids[2]);
	//fprintf(stderr, "initial: var1=%f, var2=%f, var3=%f\n", variance_cluster[0], variance_cluster[1], variance_cluster[2]);
	//fprintf(stderr, "mrf: var1=%f, var2=%f,var3=%f\n", variance_cluster[0], variance_cluster[1], variance_cluster[2]);
	//clique_pot = (double *)malloc(xdim*ydim*zdim*sizeof(double));
	//clique_pot_2 = (double *)malloc(xdim*ydim*zdim*sizeof(double));

	int max_iter = 20;
	//fprintf(stderr, "beta is:%f\n", beta);
	for(int iter=0; iter<max_iter; iter++)
	{
		//fprintf(stderr, "iteration=%d\n", iter);
		for(int k=0; k<num_clus; k++)
		{
			double mu_k = centroids[k];
			for(int ix=0;ix<xdim;ix++)          //x dim
			{
				for(int iy=0;iy<ydim;iy++)      //y dim
				{
					for(int iz=0;iz<zdim;iz++)  //z dim
					{
						double f_j  = voxvals[INDEX3D(ix, iy, iz, xdim, ydim, zdim)];
						// clear u_lm
						for(int i=0; i<num_clus*num_neib; i++)
						{
							u_lm[i] = 0.0;
						}
						for(int clus=0; clus<num_clus; clus++)
						{
							int neib_i = 0;
							for(int ia=-1; ia<2; ia++)
							{
								for(int ib=-1; ib<2; ib++)
								{
									for(int ic=-1; ic<2; ic++)
									{
										if(ix + ia < 0 || iy + ib < 0 || iz + ic < 0 || ix + ia >= xdim || iy + ib >= ydim || iz + ic >= zdim) continue;
										if(ia != 0 || ib != 0 || ic != 0)
										{
											u_lm[neib_i+clus*num_neib] = u_jk_pre[INDEX3D(ix+ia, iy+ib, iz+ic, xdim, ydim, zdim)+clus*nvox];
											neib_i++;
										}
									}// ic
								}// ib
							}// ia
							//if(INDEX3D(ix, iy, iz, xdim, ydim, zdim)%1000 == 0) fprintf(stderr, "tmp_pot1=%f\n", tmp_pot);
						}// clus

						// fuzziness power
						for(int i=0; i<num_clus*num_neib; i++)
						{
							u_lm[i]     = pow(u_lm[i], q);
							u_lm_tmp[i] = u_lm[i];
						}
						for(int i=0; i<num_neib; i++)
						{
							u_lm_tmp[i+k*num_neib] = 0;
						}
						double neib_sum = 0.0;
						for(int i=0; i<num_clus*num_neib; i++)
						{
							neib_sum += u_lm_tmp[i];
						}
						double num = pow((pow((f_j-mu_k),2.0)+beta*neib_sum),power);
						double den = 0.0;
						for(int clus = 0; clus<num_clus; clus++)
						{
							for(int i=0; i<num_clus*num_neib; i++)
							{
								u_lm_tmp[i] = u_lm[i];
							}
							for(int i=0; i<num_neib; i++)
							{
								u_lm_tmp[i+clus*num_neib] = 0;
							}
							double neib_sum = 0.0;
							for(int i=0; i<num_clus*num_neib; i++)
							{
								neib_sum += u_lm_tmp[i];
							}
							den += pow((pow((f_j-centroids[clus]),2.0)+beta*neib_sum),power);
						}
						if(isnan(den))
						{
							den = 0.00001;
						}
						if(isnan(num))
						{
							num = 0.00001;
						}
						u_jk[INDEX3D(ix, iy, iz, xdim, ydim, zdim)+k*nvox] = num/den;
						//fprintf(stderr, "u_jk[0]=%f, u_jk[1]=%f\n",u_jk[INDEX3D(ix, iy, iz, xdim, ydim, zdim)+0*nvox],u_jk[INDEX3D(ix, iy, iz, xdim, ydim, zdim)+1*nvox]);
					}// iz
				}// iy
			}// ix
		}//k class
		//fprintf(stderr, "check point1.1\n");
		for(int i=0; i<num_clus*nvox; i++)
		{
			u_jk_pre[i] = u_jk[i];
		}

		//recompute means
		for(int k=0; k<num_clus; k++)
		{
			double num = 0.0;
			double den = 0.0;
			for(int i=0; i<nvox; i++)
			{
				num += voxvals[i]*pow(u_jk[i+k*nvox],q);
				den += pow(u_jk[i+k*nvox],q);
			}
			centroids[k] = num/den;
		}
		//fprintf(stderr, "FCM: mean1=%f, mean2=%f,mean3=%f\n", centroids[0], centroids[1], centroids[2]);
	}//iter

	double pixel_class[num_clus];
	for(int i=0; i < xdim*ydim*zdim; ++i)
	{
		for(int k=0; k<num_clus; k++)
		{
			pixel_class[k] = u_jk[i+k*nvox];
		}
		labels[i] = argmax_double(pixel_class, num_clus);
	}


	free((void *)u_lm_tmp);
	free((void *)u_lm);
	free((void *)u_jk);
	free((void *)u_jk_pre);
	free_dvector((void *)centroids);
	return 0;
}

int* append_int_array(int* input_arr, int size, int num)
{
	int* output_array;
	output_array = (int *)malloc((size+1)*sizeof(int));
	if (size != 0)
	{
		for(int i=0; i<size; i++)
		{
			output_array[i] = input_arr[i];
		}
	}
	output_array[size] = num;
	free((void *)input_arr);
	return output_array;
}


int argmax_int(int* voxvals, int nvox)
{
	int max_val = voxvals[0];
	int max_i = 0;
	for(int i=0; i<nvox; i++)
	{
		if(voxvals[i]>max_val)
		{
			max_val = voxvals[i];
			max_i = i;
		}
	}
	return max_i;
}

int region_grow_cropped3Dimg(LABELTYPE* labeled_roi, int bw_ht, int bw_wid, int bw_len, double seedx, double seedy, double seedz, int num_clus)
{
	int les_reg_label;
	int wid_box    = 1;
	int count      = 0;
	int nvox       = bw_ht*bw_wid*bw_len;
	int count_temp = 0;
	int* x_ord  = NULL;
	int* y_ord  = NULL;
	int* z_ord  = NULL;
	int* temp_arr;
	int x_curr;
	int y_curr;
	int z_curr;
	LABELTYPE* flag;
	LABELTYPE* new_labels;
	flag        = (LABELTYPE *)malloc(nvox*sizeof(LABELTYPE));
	new_labels  = (LABELTYPE *)malloc(nvox*sizeof(LABELTYPE));
	temp_arr    = (int *)dvector(num_clus*sizeof(int));
	
	for(int i=0; i<num_clus; i++)
	{
		temp_arr[i] = 0;
	}
	// initialization
	for(int i=0; i < nvox; ++i)
	{
		flag[i] = 0;
		new_labels[i] = 0;
	}
	seedx = (int)seedx;
	seedy = (int)seedy;
	seedz = (int)seedz;
	
	for(int i=-wid_box; i<=wid_box; i++)
	{
		for(int j=-wid_box; j<=wid_box; j++)
		{
			for(int z=-wid_box; z<=wid_box; z++)
			{
				if(seedx+i >= bw_ht || seedy+j >= bw_wid || seedz+z >= bw_len || seedx+i < 0 || seedy+j < 0 || seedz+z < 0)// check the indexes
                	continue;
				//fprintf(stderr, "check point0\n");
				x_ord = append_int_array(x_ord, count, seedx+i);
				y_ord = append_int_array(y_ord, count, seedy+j);
				z_ord = append_int_array(z_ord, count, seedz+z);
				flag[INDEX3D((int)seedx+i, (int)seedy+j, (int)seedz+z, bw_ht, bw_wid, bw_len)] = 1;
				count ++;
				//fprintf(stderr, "check point1\n");
				for(int k=0; k<num_clus; k++)
				{
					if(labeled_roi[INDEX3D((int)seedx+i, (int)seedy+j, (int)seedz+z, bw_ht, bw_wid, bw_len)]==k)
					{
						temp_arr[k]++;
						break;
					}
				}// class k
			}// z
		}// j
	}// i
	//fprintf(stderr, "tmp_arr = %d,%d\n", temp_arr[0],temp_arr[1]);
	les_reg_label = argmax_int(temp_arr, num_clus);
	//fprintf(stderr, "label = %d\n", les_reg_label);
	//fprintf(stderr, "count = %d\n", count);
	//fprintf(stderr, "count_temp = %d\n", count_temp);
	while(count_temp<count)
	{
		//fprintf(stderr, "check point1.1\n");
		x_curr = x_ord[count_temp];
		//fprintf(stderr, "x_curr=%d\n", x_curr);
    	y_curr = y_ord[count_temp];
    	z_curr = z_ord[count_temp];
		//fprintf(stderr, "check point1.2\n");
		count_temp ++;
		//fprintf(stderr, "labeled_roi = %d\n",labeled_roi[INDEX3D(x_curr, y_curr, z_curr, bw_ht, bw_wid, bw_len)]);
		if((labeled_roi[INDEX3D(x_curr, y_curr, z_curr, bw_ht, bw_wid, bw_len)]+0.0) == (les_reg_label+0.0))
		{
			//fprintf(stderr, "in\n");
			new_labels[INDEX3D(x_curr, y_curr, z_curr, bw_ht, bw_wid, bw_len)] = 1;
			for(int ia=-1; ia<2; ia++)
			{
				for(int ib=-1; ib<2; ib++)
				{
					for(int ic=-1; ic<2; ic++)
					{
						if(x_curr+ia >=0 && y_curr+ib >= 0 && z_curr+ic >= 0 && x_curr+ia < bw_ht && y_curr+ib < bw_wid && z_curr+ic < bw_len)
						{
							if(flag[INDEX3D(x_curr+ia, y_curr+ib, z_curr+ic, bw_ht, bw_wid, bw_len)] ==0)
							{
								
								flag[INDEX3D(x_curr+ia, y_curr+ib, z_curr+ic, bw_ht, bw_wid, bw_len)] = 1;
								x_ord = append_int_array(x_ord, count, x_curr+ia);
								y_ord = append_int_array(y_ord, count, y_curr+ib);
								z_ord = append_int_array(z_ord, count, z_curr+ic);
								//fprintf(stderr, "z_ord=%d\n", z_ord[count]);
								count++;
								//fprintf(stderr, "check point3\n");
							}
						}
					}//ic
				}//ib
			}//ia
		}
	}
	for(int i=0; i < nvox; ++i)
		labeled_roi[i] = new_labels[i];
	free((void *)flag);
	free((void *)new_labels);
	free_dvector((void *)temp_arr);
	return 0;
}

int FCM_ctinfo(float* voxvals, LABELTYPE* labels, LABELTYPE* ct_label, int xdim, int ydim, int zdim, int num_clus, double beta, double q, double gamma)
{
	LABELTYPE* pixel_reg;
	int nvox     = xdim*ydim*zdim; // number of voxels
	int num_neib = 26; // number of neighborhood pixels
	double* centroids;
	double* u_jk_pre;
	double* u_jk;
	double* u_lm;
	double* u_lm_tmp;
	double* u_ct;
	double power = -1/(q-1);
	pixel_reg = (LABELTYPE *)malloc(xdim*ydim*zdim*sizeof(LABELTYPE)); // define label map
	centroids = (double *)dvector(num_clus*sizeof(double)); // class means
	u_jk_pre  = (double *)malloc(xdim*ydim*zdim*num_clus*sizeof(double));
	u_jk      = (double *)malloc(xdim*ydim*zdim*num_clus*sizeof(double));
	u_lm      = (double *)malloc(num_clus*num_neib*sizeof(double));
	u_lm_tmp  = (double *)malloc(num_clus*num_neib*sizeof(double));
	u_ct      = (double *)dvector(num_clus*sizeof(double)); // class means

	//fprintf(stderr, "class num is:%d\n", num_clus);
	k_means(voxvals, pixel_reg, xdim, ydim, zdim, num_clus, 0, 1); // k means clustering as initialization
	//for(int i=0; i < xdim*ydim*zdim; ++i)
	//	labels[i] = label_map[i];
	//return 0;
	
	for(int k=0; k<num_clus; k++)
	{
		centroids[k]     = mean_class(voxvals, nvox, pixel_reg, k);
	}//class k
	//fprintf(stderr, "initial: mean1=%f, mean2=%f, mean3=%f\n", mean_cluster[0], mean_cluster[1], mean_cluster[2]);
	//fprintf(stderr, "initial: var1=%f, var2=%f, var3=%f\n", variance_cluster[0], variance_cluster[1], variance_cluster[2]);
	//fprintf(stderr, "mrf: var1=%f, var2=%f,var3=%f\n", variance_cluster[0], variance_cluster[1], variance_cluster[2]);
	//clique_pot = (double *)malloc(xdim*ydim*zdim*sizeof(double));
	//clique_pot_2 = (double *)malloc(xdim*ydim*zdim*sizeof(double));
	//fprintf(stderr, "kmeans: u_jk[0]=%f, u_jk[1]=%f\n",u_jk[INDEX3D(ix, iy, iz, xdim, ydim, zdim)+0*nvox],u_jk[INDEX3D(ix, iy, iz, xdim, ydim, zdim)+1*nvox]);
	//fprintf(stderr, "kmeans: mean1=%f, mean2=%f,mean3=%f\n", centroids[0], centroids[1], centroids[2]);
	//fprintf(stderr, "gamma=%f\n", gamma);
	//fprintf(stderr, "q=%f\n", q);
	int max_iter = 20;
	//fprintf(stderr, "beta is:%f\n", beta);
	for(int i=0; i<xdim*ydim*zdim*num_clus; i++)
	{
		u_jk[i]     = 0;
		u_jk_pre[i] = 0;
	}
	for(int iter=0; iter<max_iter; iter++)
	{
		//fprintf(stderr, "iteration=%d\n", iter);
		for(int k=0; k<num_clus; k++)
		{
			double mu_k = centroids[k];
			for(int ix=0;ix<xdim;ix++)          //x dim
			{
				for(int iy=0;iy<ydim;iy++)      //y dim
				{
					for(int iz=0;iz<zdim;iz++)  //z dim
					{
						double f_j  = voxvals[INDEX3D(ix, iy, iz, xdim, ydim, zdim)];
						// clear u_lm
						for(int i=0; i<num_clus*num_neib; i++)
						{
							u_lm[i] = 0.0;
							u_lm_tmp[i] = 0.0;
						}
						for(int clus=0; clus<num_clus; clus++)
						{
							int neib_i = 0;
							for(int ia=-1; ia<2; ia++)
							{
								for(int ib=-1; ib<2; ib++)
								{
									for(int ic=-1; ic<2; ic++)
									{
										if(ix + ia < 0 || iy + ib < 0 || iz + ic < 0 || ix + ia >= xdim || iy + ib >= ydim || iz + ic >= zdim) continue;
										if(ia != 0 || ib != 0 || ic != 0)
										{
											u_lm[neib_i+clus*num_neib] = u_jk_pre[INDEX3D(ix+ia, iy+ib, iz+ic, xdim, ydim, zdim)+clus*nvox];
											neib_i++;
										}
									}// ic
								}// ib
							}// ia
							if(ct_label[INDEX3D(ix, iy, iz, xdim, ydim, zdim)]==clus)
							u_ct[clus] = 1;
						    else if(labels[INDEX3D(ix, iy, iz, xdim, ydim, zdim)]==num_clus-1){u_ct[num_clus-1] = 1;}//fprintf(stderr, "les=%d\n", num_clus-1);}
							else
							u_ct[clus] = 0;
						}// clus

						// smoothness
						for(int i=0; i<num_clus*num_neib; i++)
						{
							u_lm[i]     = pow(u_lm[i], q);
							u_lm_tmp[i] = u_lm[i];
						}
						for(int i=0; i<num_neib; i++)
						{
							u_lm_tmp[i+k*num_neib] = 0;
						}
						double neib_sum = 0.0;
						for(int i=0; i<num_clus*num_neib; i++)
						{
							neib_sum += u_lm_tmp[i];
						}

						//ct prior
						double ct_prior = 0.0;
						for(int c=0; c<num_clus; c++)
						{
							if(c!=k)
							{
								ct_prior += pow(u_ct[c], q);
							}
						}

						double num = pow((pow((f_j-centroids[k]),2.0)+beta*neib_sum+gamma*ct_prior),power);						
						double den = 0.0;
						for(int clus = 0; clus<num_clus; clus++)
						{
							for(int i=0; i<num_clus*num_neib; i++)
							{
								u_lm_tmp[i] = u_lm[i];
							}
							for(int i=0; i<num_neib; i++)
							{
								u_lm_tmp[i+clus*num_neib] = 0;
							}
							double neib_sum = 0.0;
							for(int i=0; i<num_clus*num_neib; i++)
							{
								neib_sum += u_lm_tmp[i];
							}

							//ct prior
							ct_prior = 0.0;
							for(int c=0; c<num_clus; c++)
							{
								if(c!=clus)
								{
									ct_prior += pow(u_ct[c], q);
								}
							}

							den += pow((pow((f_j-centroids[clus]),2.0)+beta*neib_sum+gamma*ct_prior),power);
						}
						if(isnan(den))
						{
							den = 0.00001;
						}
						if(isnan(num))
						{
							num = 0.00001;
						}
						//if(INDEX3D(ix, iy, iz, xdim, ydim, zdim)%50000 == 0) fprintf(stderr, "den=%f\n", den);
						//if(INDEX3D(ix, iy, iz, xdim, ydim, zdim)%50000 == 0) fprintf(stderr, "num=%f\n", num);
						u_jk[INDEX3D(ix, iy, iz, xdim, ydim, zdim)+k*nvox] = num/den;
					}// iz
				}// iy
			}// ix
		}//k class
		//fprintf(stderr, "check point1.1\n");
		for(int i=0; i<num_clus*nvox; i++)
		{
			u_jk_pre[i] = u_jk[i];
		}

		//recompute means
		for(int k=0; k<num_clus; k++)
		{
			double num = 0.0;
			double den = 0.0;
			for(int i=0; i<nvox; i++)
			{
				num += voxvals[i]*pow(u_jk[i+k*nvox],q);
				den += pow(u_jk[i+k*nvox],q);
				//fprintf(stderr, "FCM: den=%f\n", den);
			}
			centroids[k] = num/den;
		}
		//fprintf(stderr, "FCM: mean1=%f, mean2=%f,mean3=%f\n", centroids[0], centroids[1], centroids[2]);
		double pixel_class[num_clus];
		for(int i=0; i < xdim*ydim*zdim; ++i)
		{
			for(int k=0; k<num_clus; k++)
			{
				pixel_class[k] = u_jk[i+k*nvox];
			}
			labels[i] = argmax_double(pixel_class, num_clus);
		}
	}//iter

	


	free((void *)u_lm_tmp);
	free((void *)u_lm);
	free((void *)u_jk);
	free((void *)u_jk_pre);
	free_dvector((void *)centroids);
	return 0;
}