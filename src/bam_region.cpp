#include <algorithm>
#include "H5Cpp.h"
#include "get_reads_direct.h"
#include "region.h"
#include "bam_region.h"
#include "read.h"
#include "graph_tools.h"
#include "math_tools.h"
#include "tools.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

Bam_Region::Bam_Region()
{
	init();
}

/** constructor*/
Bam_Region::Bam_Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname)
:Region(pstart, pstop, pchr_num, pstrand, gio_fname)
{
	init();
}

Bam_Region::Bam_Region(int pstart, int pstop, int pchr_num, char pstrand)
:Region(pstart, pstop, pchr_num, pstrand)
{
	init();
}

Bam_Region::Bam_Region(int pstart, int pstop, char* pchr, char pstrand)
:Region(pstart, pstop, pchr, pstrand)
{
	init();
}

Bam_Region::Bam_Region(Region* reg)
:Region(reg)
{
	init();
}

Bam_Region::Bam_Region(Bam_Region* reg)
:Region(reg)
{
	init();
}

void Bam_Region::init()
{
	coverage = NULL;
	intron_coverage = NULL;
	reads_sorted = false;
}

/** destructor*/
Bam_Region::~Bam_Region()
{
	clear_reads();
	delete[] coverage;
	delete[] intron_coverage;
}

int Bam_Region::compute_num_paths()
{
	vector<int> initial;
	vector<int> terminal;
	initial.push_back(1);
	terminal.push_back(admat.size()-2);
	for (uint k=2; k<admat.size()-2; k++)
	{
		if (is_initial(k+1))
			initial.push_back(k);
		if (is_terminal(k+1))
			terminal.push_back(k);
	}
	int num_paths = 0;
	for (uint i=0; i<initial.size(); i++)
		for (uint k=0; k<terminal.size(); k++)
			num_paths += count_num_paths(admat, initial[i]+1, terminal[k]+1);
	return num_paths;
}

void Bam_Region::write_segment_graph(_IO_FILE*& fd)
{
	fprintf(fd, "region\t%s\t%c\t%i\t%i\n", chr, strand, start, stop);

	int num_seg = segments.size();
	fprintf(fd, "segments\t%i\n", num_seg);
	print_segments_and_coverage(fd);
	for (uint i=0; i<admat.size(); i++)
	{
		for (uint j=0; j<admat.size(); j++)
		{
			if (admat[i][j]>NEIGHBOR)
				fprintf(fd, "%i\t%i\n", i+1, j+1);
		}
	}
	fprintf(fd, "end\n");
}

int Bam_Region::compute_pair_mat()
{
	// get read pairs
	if (!reads_sorted)
	{
		sort(reads.begin(), reads.end(), CRead::compare_by_read_id);
		reads_sorted = true;
	}
	vector<int> pair_ids;
	bool no_pair_filter = true;
	
	int take_cnt = 0;
	int discard_cnt = 0; 
	int discard_strand_cnt = 0; 
	//printf("reads.size(): %i\n", (int) reads.size());
	// find consecutive reads with the same id
	for (int i=0; i<((int) reads.size())-1; i++) 
	{
		//printf("reads[%i]->read_id: %s\n", i, reads[i]->read_id);
		uint j = i+1;
		while(j<reads.size() && strcmp(reads[i]->read_id, reads[j]->read_id) == 0) 
		{
			if (((reads[i]->left && reads[j]->right) || (reads[j]->left && reads[i]->right)) && (reads[i]->reverse != reads[j]->reverse)) 
			{
				if (reads[i]->get_last_position()==-1 || reads[j]->get_last_position()==-1)
					break;
				if (false)//(reads[i]->strand[0]=='0' && reads[j]->strand[0]=='0' )
				{ 
					// discard pairs without strand information
					discard_strand_cnt++;
				}
				else if (reads[i]->get_last_position()<reads[j]->start_pos && reads[j]->start_pos-reads[i]->get_last_position()<60000) 
				{
					pair_ids.push_back(i);
					pair_ids.push_back(j);
					take_cnt++;
				}
				else if (reads[i]->start_pos>reads[j]->get_last_position() && reads[j]->get_last_position()-reads[i]->start_pos<60000)
				{
					pair_ids.push_back(i);
					pair_ids.push_back(j);
					take_cnt++;
				}
				else
				{
					if (no_pair_filter && reads[i]->start_pos<reads[j]->start_pos)
					{
						pair_ids.push_back(i);
						pair_ids.push_back(j);
						take_cnt++;
					}
					else if (no_pair_filter)
					{
						pair_ids.push_back(j);
						pair_ids.push_back(i);
						take_cnt++;
					}
					else
						discard_cnt++;
					//printf("istart:%i, ilast:%i  jstart:%i, jlast: %i\n", reads[i]->start_pos, reads[i]->get_last_position(), reads[j]->start_pos, reads[j]->get_last_position());
				}
			}
			else
				discard_cnt++;
			j++;
		}
	}

	//pair_mat = zeros(size(segments, 1));
	pair_mat.clear();
	for (uint i=0; i<segments.size(); i++)
	{
		vector<int> row(segments.size(), 0);
		pair_mat.push_back(row);
	}

	//for p = 1:size(pairs, 2)
	for (uint i=0; i<pair_ids.size(); i+=2)
	{
	//	if read_stops(pairs(1, p)+1)<read_stops(pairs(2, p)+1)
	//		r1 = pairs(1, p)+1;
	//		r2 = pairs(2, p)+1;
	//	else
	//		r1 = pairs(2, p)+1;
	//		r2 = pairs(1, p)+1;
	//	end
		int rid1 = pair_ids[i];
		int rid2 = pair_ids[i];
		int rid1_lp = reads[rid1]->get_last_position();
		int rid2_lp = reads[rid2]->get_last_position();
		if (rid1_lp>rid2_lp)
		{
			// switch reads
			int tmp = rid2;
			rid2 = rid1;
			rid1 = tmp;
			tmp = rid2_lp;
			rid2_lp = rid1_lp;
			rid1_lp = tmp;
		}
		int st1 = -1;
		int en1 = -1;
		int st2 = -1;
		int en2 = -1;
		for (uint j=0; j<segments.size(); j++)
		{
			if (segments[j].second<reads[rid1]->start_pos)
				continue;

			if (segments[j].first>rid2_lp)
				break;

			if (reads[rid1]->start_pos>=segments[j].first && reads[rid1]->start_pos<=segments[j].second)
				st1 = j;
			if (rid1_lp>=segments[j].first && rid1_lp<=segments[j].second)
				en1 = j;
			if (reads[rid2]->start_pos>=segments[j].first && reads[rid2]->start_pos<=segments[j].second)
				st2 = j;
			if (rid2_lp>=segments[j].first && rid2_lp<=segments[j].second)
				en2 = j;
		}
		//printf("r1: %i->%i, r2: %i->%i\n", reads[rid1]->start_pos, reads[rid1]->get_last_position(), reads[rid2]->start_pos, reads[rid2]->get_last_position());
		//printf("no pair_connection: %i, %i, %i, %i\n", st1, en1, st2, en2);
		if (st1>-1 && st2>-1 && en1>-1 && en2>-1)
		{
			//if (st1>=38 || st2>=38 || en1>=38 || en2>=38)
			//	printf("add pair_connection: %i, %i, %i, %i\n", st1, en1, st2, en2);
			pair_mat[st1][st2]++;
			pair_mat[st2][st1]++;

			pair_mat[en1][en2]++;
			pair_mat[en2][en1]++;

			pair_mat[st1][en2]++;
			pair_mat[en2][st1]++;

			pair_mat[en1][st2]++;
			pair_mat[st2][en1]++;
		}
		//else if (st1>=38 || st2>=38 || en1>=38 || en2>=38)
		//	printf("no pair_connection: %i, %i, %i, %i\n", st1, en1, st2, en2);
	//	% find the segments where the left read starts and right read ends
	//	% and increase pair counter between these segments
	//	for x = 1:4
	//		if x==1 % start <-> stop
	//			s1 = find(segments(:,1)<=read_starts(r1)&segments(:,2)>=read_starts(r1));
	//			s2 = find(segments(:,1)<=read_stops(r2)&segments(:,2)>=read_stops(r2));
	//		elseif x==2 % start <-> start
	//			s1 = find(segments(:,1)<=read_starts(r1)&segments(:,2)>=read_starts(r1));
	//			s2 = find(segments(:,1)<=read_starts(r2)&segments(:,2)>=read_starts(r2));
	//		elseif x==3 % stop <-> stop
	//			s1 = find(segments(:,1)<=read_stops(r1)&segments(:,2)>=read_stops(r1));
	//			s2 = find(segments(:,1)<=read_stops(r2)&segments(:,2)>=read_stops(r2));
	//		elseif x==4 % stop <-> start
	//			s1 = find(segments(:,1)<=read_stops(r1)&segments(:,2)>=read_stops(r1));
	//			s2 = find(segments(:,1)<=read_starts(r2)&segments(:,2)>=read_starts(r2));
	//		end

	//		if isempty(s1) || isempty(s2)
	//			continue
	//		end

	//		assert(length(s1)==1)
	//		assert(length(s2)==1)
	//		
	//		pair_mat(s1, s2) = pair_mat(s1, s2)+1;
	//		if s1~=s2
	//			pair_mat(s2, s1) = pair_mat(s2, s1)+1;
	//		end
	//	end
	//end
	}
	return 0;
}

void Bam_Region::clear_reads()
{
	//for (int i=0; i<all_reads.size(); i++)
	//	delete all_reads[i];
	all_reads.clear();
	// reads is a subset of all_reads, therefore the 
	// destructor for each read has already been called
	// only the pointers have to be removed
	reads.clear(); 
}


void Bam_Region::get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter)
{
	get_reads(bam_files, num_bam_files, intron_len_filter, filter_mismatch, exon_len_filter, false);
}

void Bam_Region::get_reads(char** bam_files, int num_bam_files, int intron_len_filter, int filter_mismatch, int exon_len_filter, bool mm_filter)
{

	char* reg_str = get_region_str();
	int subsample = 1000;
	for (int i=0; i<num_bam_files; i++)
	{
    	char* fn_bam = bam_files[i];
	    fprintf(fd_out, "getting reads from file %i: %s\n", i, fn_bam);
		get_reads_from_bam(fn_bam, reg_str, &all_reads, strand, subsample);
		fprintf(fd_out, "#reads: %lu\n", all_reads.size());
	}
	delete[] reg_str; 
	fprintf(fd_out, "number of reads (not filtered): %lu\n", all_reads.size());
	/* filter reads
	* **************/
	//int exon_len_filter = 8;
	//int filter_mismatch = 1;
	//int intron_len_filter = 100000;
	int intron_filter = 0;
	int exon_filter = 0;
	int mismatch_filter = 0;
	int multi_filter= 0;
	for (uint i=0; i<all_reads.size(); i++)
	{
		bool take = true;
		if (all_reads[i].max_intron_len()>=intron_len_filter)
		{
			intron_filter++;
			take = false;
		}
		if (all_reads[i].min_exon_len()<=exon_len_filter)
		{
			exon_filter++;
			take = false;
		}
		if (all_reads[i].get_mismatches()>filter_mismatch)
		{
			mismatch_filter++;
			take = false;
		}
		if (all_reads[i].multiple_alignment_index!=0 && mm_filter)
		{
			multi_filter++;
			take = false;
		}

	    if (take)
	        reads.push_back(&all_reads[i]);
	}
	fprintf(fd_out, "number of reads: %lu\n", reads.size());
	fprintf(fd_out, "Filter: intron_len:%i, min_exon_len:%i, mismatch:%i, multi_mapper:%i\n", intron_filter, exon_filter, mismatch_filter, multi_filter);
}



void Bam_Region::compute_coverage() 
{
	int num_pos = stop-start+1;
	if (!coverage)
		coverage = new uint32_t[num_pos];

	memset(coverage, 0, num_pos*sizeof(uint32_t));
	//for (int i=0; i<num_pos; i++)
	//	coverage[i] = 0;
	for (uint i=0; i<reads.size(); i++)
	{
		// add read contribution to coverage
		reads[i]->get_coverage(start, stop, coverage);
	}
}

float Bam_Region::get_coverage_global(int pstart, int pstop)
{
	if (!coverage)
		compute_coverage();
	
	assert(pstart<=pstop);

	int sum = 0;
	for (int i = pstart; i<=pstop; i++)
	{
		if (!(i>=start && i<=stop))
		{
			fprintf(stdout, "start%i pstart%i stop%i pstop%i\n",start, pstart, stop, pstop);	
		}
		assert(i>=start && i<=stop);
		sum+=coverage[i-start];
	}
	return ((float) sum)/(pstop-pstart+1);
}

float Bam_Region::get_coverage_seg(int i)
{
	assert(i<(int) segments.size());
	if (seg_cov.size()==0)
	{
		compute_seg_cov();
	}
	return seg_cov[i];
}

void Bam_Region::compute_seg_cov()
{
	seg_cov.clear();
	for (uint s=0; s<segments.size(); s++)
		seg_cov.push_back(get_coverage_global(segments[s].first, segments[s].second));
}

void Bam_Region::compute_intron_coverage() 
{
	int num_pos = stop-start+1;
	if (!intron_coverage)
		intron_coverage = new uint32_t[num_pos];

	for (int i=0; i<num_pos; i++)
		intron_coverage[i] = 0;
	for (uint i=1; i<unique_introns.size(); i++)
	{
		int from_pos = std::max(start, unique_introns[i].first);
        int to_pos = std::min(stop, unique_introns[i].second);
		for (int j=from_pos; j<to_pos; j++)
		{
			intron_coverage[j] += intron_counts[i];
		}
	}
}

int Bam_Region::get_intron_conf(int intron_start, int intron_stop)
{
	//printf("get_intron_conf for %i -> %i\n", intron_start, intron_stop);
	for (uint i=1; i<unique_introns.size(); i++)
	{
		//if (unique_introns[i].first==intron_start || unique_introns[i].second==intron_stop)
		//printf("found intron %i -> %i\n", unique_introns[i].first, unique_introns[i].second);
		if (unique_introns[i].first==intron_start && unique_introns[i].second==intron_stop)
			return intron_counts[i];
	}
	return 0;
}

int Bam_Region::get_read_starts(int from, int to)
{
	int ret = 0;
	for (uint i=0; i<reads.size(); i++)
		if (reads[i]->start_pos>=from && reads[i]->start_pos<=to)
			ret++;
	return ret;
}
int Bam_Region::get_read_ends(int from, int to)
{
	int ret = 0;
	for (uint i=0; i<reads.size(); i++)
		if (reads[i]->get_last_position()>=from && reads[i]->get_last_position()<=to)
			ret++;
	return ret;
}
void Bam_Region::find_tss_and_cleave(vector<int>* pos, vector<int>* starts, vector<int>* stops, float pval)
{
	// process read starts and ends
	vector<int> rstarts;
	vector<int> rends;
	for (uint i=0; i<reads.size(); i++)
	{
		rstarts.push_back(reads[i]->start_pos);
		rends.push_back(reads[i]->get_last_position());
	}
	sort(rstarts.begin(), rstarts.end());
	sort(rends.begin(), rends.end());

	if (rstarts.size()==0)
		return;

	int* rs = &rstarts[0];
	int* re = &rends[0];

	int num_pos = pos->size();
	for (int i=0; i<num_pos-1; i++)
	{
		if (pos->at(i)==pos->at(i+1))
			continue;
		
		int win=30;
		vector<int> bins;
		vector<int> start_cnts;
		vector<int> end_cnts;
		for (int j=pos->at(i); j<pos->at(i+1)-win; j+=win)
		{
			bins.push_back(j);
			start_cnts.push_back(0);
			end_cnts.push_back(0);
			while (rs<&rstarts.back() && *rs<j+win)
			{
				rs++;
				if (*rs>=j)
					start_cnts.back()++;
			}
			//printf("%i, %i\n", start_cnts.back(), get_read_starts(j+1, j+win));
			while (re<&rends.back() && *re<j+win)
			{
				re++;
				if (*re>=j)
					end_cnts.back()++;
			}
			//printf("%i, %i\n", end_cnts.back(), get_read_ends(j+1, j+win));
		}
		for (int b=0; b<((int) bins.size())-1; b++)
		{
			//int cnt1 = get_read_starts(j, j+win-1);
			//int cnt2 = get_read_starts(j+win, j+2*win-1);

			//int cnt3 = get_read_ends(j, j+win-1);
			//int cnt4 = get_read_ends(j+win, j+2*win-1);

			int cnt1 = start_cnts[b];
			int cnt2 = start_cnts[b+1];

			int cnt3 = end_cnts[b];
			int cnt4 = end_cnts[b+1];

			float p_val1 = bino_cfd(cnt1, cnt1+cnt2, 0.5);
			float p_val2 = bino_cfd(cnt4, cnt3+cnt4, 0.5);
			if (p_val1<pval)
			{
				//printf("add tss%i: cnt1:%i, cnt2:%i, pval:%.6f\n", j+win, cnt1, cnt2, p_val1);
				pos->push_back(bins[b]+win);
				starts->push_back(bins[b]+win);
			}
			else if (p_val2<pval)
			{
				//printf("add cleave%i: cnt1:%i, cnt2:%i, pval:%.6f\n", j+win, cnt1, cnt2, p_val2);
				pos->push_back(bins[b]+win);
				stops->push_back(bins[b]+win);
			}
		}
	}
}


void Bam_Region::generate_segment_graph(float seg_filter, float tss_pval)
{
	compute_intron_list();
	vector<int> pos;
	vector<int> starts;
	vector<int> stops;

	// add splice sites, tss, and tts from transcripts
	for (uint i=0; i<transcripts.size(); i++)
	{
		starts.push_back(transcripts[i].front().first-1);
		stops.push_back(transcripts[i].back().second);
		for (uint j=0; j<transcripts[i].size(); j++)
		{
			pos.push_back(transcripts[i][j].first-1);
			pos.push_back(transcripts[i][j].second);
		}
	}
	//printf("pos from anno:");
	//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]);printf("\n");

	if (reads.size()==0 && pos.size()==0)
	{
		segment seg(start, stop);
		segments.push_back(seg);
	}
	else
	{
		pos.push_back(start-1);
		for (uint i=0; i<unique_introns.size(); i++)
		{
			if (unique_introns[i].first>start && unique_introns[i].first<stop)
				pos.push_back(unique_introns[i].first-1);// last exonic position
			if (unique_introns[i].second>start && unique_introns[i].second<stop)
				pos.push_back(unique_introns[i].second);// the following segment starts with the first exonic position
		}

		//printf("pos from spliced reads: "); 
		//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]); printf("\n");

		// add segment boundaries where coverage drops to zero
		compute_coverage();
		int num_pos = stop-start+1;
		for (int i=0; i<num_pos-1; i++)
		{
			if (coverage[i]==0 && coverage[i+1]>0)
				pos.push_back(i+1+start-1);// start next segment with the first exonic position
			else if (coverage[i]>0 && coverage[i+1]==0)
				pos.push_back(i+start);// this is the last exonic position
		}

		sort(pos.begin(), pos.end());
		pos.push_back(stop);
		
		//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]);
		//printf("\n");

		find_tss_and_cleave(&pos, &starts, &stops, tss_pval);

		sort(pos.begin(), pos.end());

		//for (int i=0; i<pos.size(); i++) printf("%i, ", pos[i]);
		//printf("\n");
		
		//for (int i=0; i<starts.size(); i++) printf("%i, ", starts[i]);
		//printf("\n");
		//for (int i=0; i<stops.size(); i++) printf("%i, ", stops[i]);
		//printf("\n");

		
		// create segments
		for (int i=0; i<((int)pos.size())-1; i++)
		{
			if (pos[i]==pos[i+1])
				continue;
			segment seg(pos[i]+1, pos[i+1]);
			segments.push_back(seg);
			//printf("segments: %i->%i\n", pos[i]+1, pos[i+1]);
		}

		compute_pair_mat();
		compute_admat(starts, stops);
		// filter segments without coverage
		vector<segment> seg_filtered;
		for (uint i=0; i<segments.size(); i++)
		{
			if (!(segments[i].first<=segments[i].second))
			{
				
				print_segments_and_coverage(stdout); 
				assert(false);
			}
			float cov = 0;
			for (int j=segments[i].first; j<=segments[i].second; j++)
				if (coverage[j-start]>0)
					cov++;

			cov = cov/(segments[i].second-segments[i].first+1);
			//float cov = get_coverage_global(segments[i].first, segments[i].second);
			//printf("%i->%i  %.4f, %.4f\n", segments[i].first, segments[i].second, cov, seg_filter);
			
			// check if segment is part of any transcript
			bool keep = cov>seg_filter || is_annotated(i);


			// check if removal of this segment cuts a path between to 
			// paired end reads

			if (!keep && segments[i].second-segments[i].first<1000)
			{

				// remove all connections to the current segment
				vector<vector<float> > admat_wo = admat;
				for (uint j=0; j<admat.size(); j++)
				{
					admat_wo[i+1][j] = NO_CONNECTION;
					admat_wo[j][i+1] = NO_CONNECTION;
				}


				for (uint j=0; j<i; j++)
				{
					for (uint k=i+1; k<pair_mat.size(); k++)
					{
						if (j==k)
							continue;
						if (pair_mat[j][k]>0)
						{
							//printf("seg: %i checking segments: %i->%i\n", i, j, k);
							int np2 = count_num_paths(admat_wo, j+1, k+1);
							if (np2>0)
								continue;
							int np1 = count_num_paths(admat, j+1, k+1);
							if (np2==0 && np1>0)
							{
								// had only one connection via the current segment
								//printf("pairs: keep segment: %i->%i\n", segments[i].first, segments[i].second);
								keep = true;
							}
						}
					}
				}
			}


			if (keep) 
			{
				seg_filtered.push_back(segments[i]);	
			}
			//else
			//	printf("discard segment[%i]: %i->%i\n", i, segments[i].first, segments[i].second);
		}

		//printf("\t%i (%i) segments passed filter\n", (int)seg_filtered.size(), (int) segments.size());
		if (seg_filtered.size()>0)
		{
			segments.clear();
			segments = seg_filtered;
		}
		compute_admat(starts, stops);
		// filter segments without coverage

		seg_filtered.clear();
		for (uint i=0; i<segments.size(); i++)
		{

			bool no_parents = get_parents(i+1, false).empty();
			bool no_children = get_children(i+1, false).empty();	
			//printf("segment: %i->%i, no_parents, %i, no_children: %i\n", segments[i].first, segments[i].second, no_parents, no_children);

			if (is_annotated(i) || !no_parents || !no_children || segments[i].second-segments[i].first>150)
			{
				// discard segment even if well covered
				seg_filtered.push_back(segments[i]);
				//printf("no parents, no children: discard\n");
			}
		}
		if (seg_filtered.size()>0)
		{
			segments.clear();
			segments = seg_filtered;
		}
	}
	compute_admat(starts, stops);


	// remove short tss and tts segment next to splice sites
	// those originate most likely from wrong alignments of 
	// reads that should be spliced
	uint seg=0;
	while (seg<segments.size())
	{
		//printf("checking segment: %i->%i\n", segments[seg].first, segments[seg].second);
		bool discard = false;
		if (is_annotated(seg))
		{
			//printf("is annotated\n");
			seg++;
			continue;
		}

		if (segments[seg].second-segments[seg].first>20)
		{
			//printf("is long\n");
			seg++;
			continue;
		}
		if (is_initial(seg+1))
		{
			//printf("is initial, next is acceptor: %i\n", is_acceptor_ss(seg+1+1));
			if (seg+1<segments.size() && segments[seg+1].first==segments[seg].second+1 && is_acceptor_ss(seg+1+1))// +1 because admat is 1 based
			{
				//printf("discard\n");
				discard = true;
			}
		}
		else if(is_terminal(seg+1))
		{
			//printf("is terminal, previous is donor: %i\n", (int) is_donor_ss(seg-1+1));
			if (seg>0 && segments[seg-1].second+1==segments[seg].first && is_donor_ss(seg-1+1))// +1 because admat is 1 based 
			{
				//printf("discard\n");
				discard = true;
			}
		}

		if (discard)
		{
			vector<segment> seg_filtered;
			for (uint j=0; j<segments.size(); j++)
			{
				if (j!=seg)
					seg_filtered.push_back(segments[j]);
			}
			segments = seg_filtered;
			compute_admat(starts, stops);
		}
		else
			seg++;
	}

	// compute transcript paths
	for (uint i=0; i<transcripts.size(); i++)
	{
		vector<int> path; 
		uint seg = 0; 
		while (segments[seg].first!=transcripts[i][0].first)
		{
			seg++;
			assert(seg<segments.size());
		}
		assert(seg<segments.size());

		path.push_back(seg);
		uint exon = 0;
		while (exon<transcripts[i].size())
		{
			assert(segments[seg].first==transcripts[i][exon].first);
			while (segments[seg].second != transcripts[i][exon].second)
			{
				assert(seg+1<segments.size());
				// make sure segments are connected
				assert(segments[seg].second+1==segments[seg+1].first);
				seg++;
				//add segment to the path
				path.push_back(seg);
			}

			exon++;
			if (exon<transcripts[i].size())
			{
				// find segment
				while (segments[seg].first!=transcripts[i][exon].first)
				{
					seg++; 
					assert(seg<segments.size());
				}
				assert(segments[seg].first==transcripts[i][exon].first);
				path.push_back(seg);
				// assert that intron connection exists
			}
		}
		transcript_paths.push_back(path);
		// print
		//for (int exon=0; exon<transcripts[i].size(); exon++)
		//	printf("exon: %i->%i\n", transcripts[i][exon].first, transcripts[i][exon].second );

		//for (int seg=0; seg<segments.size(); seg++)
		//	printf("segment%i: %i->%i\n", seg, segments[seg].first, segments[seg].second);

		//printf("path:  ");
		for (uint p=0; p<path.size(); p++)
		{
			//printf("%i\n", path[p]);
			if (p>0)
			{
				//printf("\t%.3f\n", admat[path[p-1]+1][path[p]+1]);
				assert(path[p-1]<path[p]);
				assert(admat[path[p-1]+1][path[p]+1]>=NEIGHBOR);
			}
		}
		//printf("\n");
		

		// recompute pair mat
		compute_pair_mat();
	}
}
bool Bam_Region::is_annotated(int i)
{
	bool ret = false;
	for (uint k=0; k<transcripts.size(); k++)
		for (uint j=0; j<transcripts[k].size(); j++)
			if (transcripts[k][j].first<=segments[i].first && transcripts[k][j].second>=segments[i].second)
				ret = true;
	return ret;
}
bool Bam_Region::is_donor_ss(int i)
{
	bool ret = false;
	for (uint j=i+1; j<=segments.size(); j++)
	{
		if (admat[i][j]>NO_CONNECTION && segments[i-1].second<segments[j-1].first-1)
			ret = true;
	}
	return ret;
}
bool Bam_Region::is_acceptor_ss(int i)
{
	bool ret = false;
	for (int j=1; j<i; j++)
	{
		if (admat[j][i]>NO_CONNECTION && segments[j-1].second<segments[i-1].first-1)
			ret = true;
	}
	return ret;
}
bool Bam_Region::is_initial(int i)
{
	assert(i < (int) admat[0].size());
	return admat[0][i]>NO_CONNECTION || i==1;
}
bool Bam_Region::is_terminal(int i)
{
	int num_seg = segments.size();
	assert(i < (int) admat[num_seg].size());
	return admat[num_seg+1][i]>NO_CONNECTION || i==num_seg;
}
vector<int> Bam_Region::find_max_path()
{
	// find path with highest expression
	int max_seg = -1;
	float max_cov = -1;
	for (uint i=0; i<segments.size(); i++)
	{
		float cov = get_coverage_global(segments[i].first, segments[i].second);
		if (cov>max_cov)
		{
			max_cov=cov;
			max_seg=i;
		}
	}
	vector<int> path;
	path.push_back(max_seg);


	// search upstream
	int i=max_seg;
	while (!is_initial(i+1))
	{
		vector<int> parents = get_parents(i, false);	
		assert(parents.size()>0);
		float max_cov = -3;
		int next = -1;
		for (uint j=0; j<parents.size(); j++)
		{
			int k = parents[j];
			float cov=-1;
			if (admat[i][k]==NEIGHBOR)
			{
				float ci = get_coverage_seg(i);
				float ck = get_coverage_seg(k);
				cov = std::min(ci, ck);
			}
			else if (admat[i][k]==NO_CONNECTION)
			{
				fprintf(stderr, "This is a bug: no connection from %i->%i\n", i, k);
			}
			else
			{
				cov = admat[i][k];
			}
			if (cov>max_cov)
			{
				max_cov = cov;
				next = j;
			}
		}
		assert(next>-1);
		i = next;
		fprintf(stdout, "%i  ", next);
		path.push_back(next);
	}
	// search downstream
	i=max_seg;
	while (!is_terminal(i))
	{
		vector<int> children = get_children(i, false);	
		float max_cov = -1;
		int next = -1;
		for (uint j=0; j<children.size(); j++)
		{
			int k = children[j];
			float cov=-1;
			if (admat[i][k]==NEIGHBOR)
			{
				float ci = get_coverage_seg(i);
				float ck = get_coverage_seg(k);
				cov = std::min(ci, ck);
			}
			else
			{
				cov = admat[i][k];
			}
			if (cov>max_cov)
			{
				max_cov = cov;
				next = j;
			}
		}
		i = next;
		fprintf(stdout, "%i  ", next);
		path.push_back(next);
	}
	// sort path
	sort(path.begin(), path.end());

	return path;
}
void Bam_Region::add_bias_counts(vector<int>* vec)
{
	//int num_bins = vec->size();
	vector<float> exonic_cov;

	vector<int> path = find_max_path();
	//printf("\n");
	for (uint i=0; i<path.size(); i++)
	{
		int j=path[i];
		for (int p=segments[j].first; p<segments[j].second; p++)
		{
			//float cov = get_coverage_global(p, p);
			exonic_cov.push_back(p);
	//		printf("%.1f ", cov);
		}
	}
	//printf("\n");

	//int exonic_len = exonic_cov.size();
	//float step = exonic_len/num_bins;

}

void Bam_Region::print_segments_and_coverage(_IO_FILE*& fd)
{
	for (uint i=0; i<segments.size(); i++) 
	{
		float cov = get_coverage_global(segments[i].first, segments[i].second);
		fprintf(fd, "%i\t%i\t%.4f\n", segments[i].first, segments[i].second, cov);
	}
}

void Bam_Region::init_admat(int num_seg)
{
	admat.resize(num_seg+2); // source and sink node
	for (uint i=0; i<admat.size(); i++)
	{
		admat[i].resize(num_seg+2);
		for (uint j=0; j<admat.size(); j++)
			admat[i][j] = NO_CONNECTION;
	}
}

void Bam_Region::compute_admat(vector<int> starts, vector<int> stops)
{
	// check that segments are sorted
	for (uint i=1; i<segments.size(); i++)
		assert(segments[i].first>=segments[i-1].second);

	// allocate memory and set all entries to NO_CONNECTION
	int num_seg = segments.size();
	init_admat(num_seg);

	// add connections from RNA-seq introns
	update_coverage_information();

	// connect neighboring segments
	for (uint j=1; j<segments.size(); j++)
	{
		if (segments[j-1].second+1==segments[j].first)
		{
			admat[j][j+1] = NEIGHBOR;
			admat[j+1][j] = NEIGHBOR;
			//printf("connecting neighbor (%i->%i) %i, %i\n", segments[j-1].second, segments[j].first, j, j+1);
		}
		//else
			//printf("not connecting neighbor (%i->%i) %i, %i\n", segments[j-1].second, segments[j].first, j, j+1);
	}
	// connect nodes to source and sink
	for (uint j=0; j<segments.size(); j++)
	{
		vector<int> parents = get_parents(j+1, false);	
		vector<int> children = get_children(j+1, false);	
		bool is_start=false;
		for (uint k=0; k<starts.size(); k++)
			if(segments[j].first==starts[k]+1)
				is_start = true;
		bool is_stop=false;
		for (uint k=0; k<stops.size(); k++)
			if(segments[j].second==stops[k])
				is_stop = true;

		//printf("%i->%i, start%i, stop%i parents:%i, children:%i\n", segments[j].first, segments[j].second, is_start, is_stop, (int)parents.size(), (int)children.size());

		if (parents.size()==0 || is_start)
		{
			admat[0][j+1] = CONNECTION;
			admat[j+1][0] = CONNECTION;
		}
		if (children.size()==0 || is_stop)
		{
			admat[num_seg+1][j+1] = CONNECTION;
			admat[j+1][num_seg+1] = CONNECTION;
		}
	}
}
void Bam_Region::update_coverage_information()
{
	compute_intron_list();

	// add connections from transcripts
	vector<segment> trans_introns;
	for (uint i=0; i<transcripts.size(); i++)
	{
		int num_exons = transcripts[i].size();
		for (int j=0; j<num_exons-1; j++)
		{
			segment* seg = new segment();
			// intron coordinates should be: 
			//       |-------|
			//XXXXXXX---------XXXXXXX
			// exon1   intron  exon2
			seg->first = transcripts[i][j].second+1;
			seg->second = transcripts[i][j+1].first-1;
			trans_introns.push_back(*seg);
			delete seg;
		}
	}
	// match transcript introns with RNA-seq introns
	int match_cnt = 0;
	for (uint i=0; i<trans_introns.size(); i++)
	{
		bool matched = false;
		for (uint j=0; j<unique_introns.size(); j++)
		{
			if (trans_introns[i].first==unique_introns[j].first && trans_introns[i].second==unique_introns[j].second)
			{
				matched = true;
				break;
			}
		}
		if (!matched)
		{
			unique_introns.push_back(trans_introns[i]);
			intron_counts.push_back(0);
		}
		else
			match_cnt++;
	}
	if (trans_introns.size()>5 && unique_introns.size()>trans_introns.size()+5 && match_cnt<2)
	{
		printf("Warning: annotated introns do not match RNA-seq introns\n");
		for (uint i=0; i<unique_introns.size(); i++)
		{
			printf("%i->%i (%i)\n", unique_introns[i].first, unique_introns[i].second, intron_counts[i]);
		}
	}

	// remove previous coverage information
	for (uint i=0; i<admat.size(); i++)
		for (uint j=0; j<admat.size(); j++)
			if (admat[i][j]>=CONNECTION)
				admat[i][j] = CONNECTION;

	// add counts for introns
	for (uint i=0; i<unique_introns.size(); i++)
	{
		int start = unique_introns[i].first-1;
		int stop = unique_introns[i].second+1;
		int seg1=-1;
		int seg2=-1;
		for (uint j=0; j<segments.size(); j++)
		{
			if (segments[j].second==start)
			{
				seg1 = j+1;// one based because of source node
			}
			if (segments[j].first==stop)
			{
				seg2 = j+1;// one based because of source node
			}
		}
		//assert(seg1>=0 && seg2>=0);
		if (!(seg1>=0 && seg2>=0))
			continue;
		//{
		//	printf("intron missing: %i->%i\n", start, stop);
		//}
		//printf("intron: %i->%i %i->%i (%i)\n", start, stop, seg1, seg2, intron_counts[i]);
		if (intron_counts[i]==0)
		{
			assert(admat[seg1][seg2]<=0);
			assert(admat[seg2][seg1]<=0);
		}
		admat[seg2][seg1] = intron_counts[i];
		admat[seg1][seg2] = intron_counts[i];
	}
}
vector<int> Bam_Region::get_parents(int node, bool no_neighbors)
{
	int thresh;
	if (no_neighbors)
		thresh = CONNECTION;
	else
		thresh = NEIGHBOR;

	vector<int> parents;
	for (int j=1; j<node; j++)
	{
		if (admat[j][node]>=thresh)
			parents.push_back(j);
	}
	return parents;
}
vector<int> Bam_Region::get_children(int node, bool no_neighbors)
{
	int thresh;
	if (no_neighbors)
		thresh = CONNECTION;
	else
		thresh = NEIGHBOR;

	vector<int> children;
	for (uint j=node+1; j<=segments.size(); j++)
	{
		if (admat[j][node]>=thresh)
			children.push_back(j);
	}
	return children;
}
int Bam_Region::read_HDF5(char* filename, int graph_idx)
{
	typedef struct meta_info {
	int    num_graphs;
	} meta_info;
	
	try
	{
		// suppress printing of exceptions to allow for correct handling
		Exception::dontPrint();

		const H5std_string FILE_NAME(filename);
		H5File* file;
		try 
		{
			// this failes if file does not exist
			file = new H5File(FILE_NAME, H5F_ACC_RDWR);// read and write access to existing file
		}
		catch (Exception error)
		{
			fprintf(stderr, "Could not open file: %s\n", filename); 
			return -1;
		}

		// read the number of graphs currently stored in the file from the file's 
		// meta data 
		CompType minfo(sizeof(meta_info));
		const H5std_string member_name( "num_graphs" );
		minfo.insertMember( member_name, HOFFSET(meta_info, num_graphs), PredType::NATIVE_INT);
		
		meta_info meta[1];
		meta[0].num_graphs = -1;
		DataSet*	metadata;

		const	H5std_string DATASET_NAME( "Graph_meta_info" );
		
		try 
		{  // to determine if the dataset exists in the file
			metadata = new DataSet(file->openDataSet(DATASET_NAME));
			metadata->read( meta, minfo );
		}
		catch( FileIException not_found_error )
		{
			fprintf(stderr, "Dataset not found\n");
			return -1;
		}
		
		//printf("number of graphs in file: %i\n", meta[0].num_graphs);

		if (meta[0].num_graphs<graph_idx)
		{
			fprintf(stderr, "index %i exceeds number of graphs %i\n", graph_idx, meta[0].num_graphs);
			return -1;
		}

		char group_name[1000];
		sprintf(group_name, "/Graph_%i", graph_idx);
		Group* group;
		try
		{
			group = new Group(file->openGroup(group_name));
		}
		catch (GroupIException error)
		{
			fprintf(stderr, "group %s not found\n", group_name);
			return -1;
		}

		// Create string datatype
    	StrType tid1(0, H5T_VARIABLE);
    	if(H5T_STRING!=H5Tget_class(tid1.getId()) || !H5Tis_variable_str(tid1.getId()))
       		printf("this is not a variable length string type!!!");
	

		{
			// read chromosome name and strand
			char d_name[1000];
			sprintf(d_name, "%s/region_str", group_name);
			DataSet dataset = file->openDataSet(d_name);

    		DataType dtype = dataset.getDataType();
    		assert(H5Tequal(H5Tget_native_type(dtype.getId(), H5T_DIR_DEFAULT), tid1.getId()));

   			DataSpace sid1 = dataset.getSpace();
			int rank = sid1.getSimpleExtentNdims();
			assert(rank==1);
			hsize_t		dims1[rank];
			sid1.getSimpleExtentDims( dims1, NULL);
			
			char* reg_str[dims1[0]];


			dataset.read((void*) reg_str, dtype);

			assert(strlen(reg_str[0])>0&strlen(reg_str[0])<1000);
			assert(strlen(reg_str[1])>0);
			strand = reg_str[1][0];
			chr = new char[strlen(reg_str[0])+1];
			for (int i=0; i<strlen(reg_str[0]); i++)
			{
				if (reg_str[0][i]=='-' || reg_str[0][i]==':')
					reg_str[0][i]= '\t';
			}
			int num_read = sscanf(reg_str[0], "%s\t%i\t%i", chr, &start, &stop);
			if (num_read!=3)
			{
				fprintf(stderr, "bam_region: Error parsing line: %s (num_read:%i), chr:%s\n", reg_str[0], num_read, chr);
				return -1;
			}
			free(reg_str[0]);
			free(reg_str[1]);

    		/* Close Dataset */
    		dataset.close();
		}

		try
		{
			// read list of segments from file
			//
			segments.clear();

			char d_name[1000];
			sprintf(d_name, "%s/segments", group_name);
			DataSet* dataset_seg = new DataSet(file->openDataSet(d_name));

			DataType dtype_seg = dataset_seg->getDataType();
			assert(H5Tequal(H5Tget_native_type(dtype_seg.getId(), H5T_DIR_DEFAULT), PredType::NATIVE_INT.getId()));

			DataSpace dataspace_seg = dataset_seg->getSpace();
			assert(dataspace_seg.getSimpleExtentNdims()==2);

			hsize_t dims_seg[2];
			dataspace_seg.getSimpleExtentDims(dims_seg, NULL);
			assert(dims_seg[1]==3);

			int seg[dims_seg[0]][3];
			dataset_seg->read(seg, PredType::NATIVE_INT);

			for (int i=0; i<dims_seg[0]; i++)
			{
				segment tmp(seg[i][0], seg[i][1], seg[i][2]);
				segments.push_back(tmp);
			}
			delete dataset_seg;
		}
		catch( FileIException not_found_error )
		{
			fprintf(stderr, "segments not found\n");
		}

		if (segments.size()>0)
		{
			try 
			{  // to determine if the dataset exists in the file
				char d_name[1000];
				sprintf(d_name, "%s/admat_idx1", group_name);
				DataSet* dataset_idx1 = new DataSet(file->openDataSet(d_name));

				sprintf(d_name, "%s/admat_idx2", group_name);
				DataSet* dataset_idx2 = new DataSet(file->openDataSet(d_name));

				sprintf(d_name, "%s/admat_val", group_name);
				DataSet* dataset_val = new DataSet(file->openDataSet(d_name));

				// check data types
	    		DataType dtype_idx1 = dataset_idx1->getDataType();
	    		DataType dtype_idx2 = dataset_idx2->getDataType();
	    		DataType dtype_val = dataset_val->getDataType();
    		
    			assert(H5Tequal(H5Tget_native_type(dtype_idx1.getId(), H5T_DIR_DEFAULT), PredType::NATIVE_INT.getId()));
    			assert(H5Tequal(H5Tget_native_type(dtype_idx2.getId(), H5T_DIR_DEFAULT), PredType::NATIVE_INT.getId()));
    			assert(H5Tequal(H5Tget_native_type(dtype_val.getId(), H5T_DIR_DEFAULT), PredType::NATIVE_FLOAT.getId()));

				// get data space and assert dimensions
   				DataSpace dataspace_idx1 = dataset_idx1->getSpace();
   				DataSpace dataspace_idx2 = dataset_idx2->getSpace();
   				DataSpace dataspace_val = dataset_val->getSpace();

				assert(dataspace_idx1.getSimpleExtentNdims()==1);
				assert(dataspace_idx2.getSimpleExtentNdims()==1);
				assert(dataspace_val.getSimpleExtentNdims()==1);

				hsize_t dims_idx1[1];
				dataspace_idx1.getSimpleExtentDims(dims_idx1, NULL);
				
				hsize_t dims_idx2[1];
				dataspace_idx2.getSimpleExtentDims(dims_idx2, NULL);

				hsize_t dims_val[1];
				dataspace_val.getSimpleExtentDims(dims_val, NULL);
				
				assert(dims_idx1[0]==dims_idx2[0]);
				assert(dims_idx1[0]==dims_val[0]);

				// init buffer
				int idx1[dims_idx1[0]];
				int idx2[dims_idx1[0]];
				int val[dims_idx1[0]];

				// read from disk
				dataset_idx1->read(idx1, PredType::NATIVE_INT);
				dataset_idx2->read(idx2, PredType::NATIVE_INT);
				dataset_val->read(val, PredType::NATIVE_FLOAT);

				init_admat(segments.size());

				for (int i=0; i<dims_idx1[0]; i++)
				{
					assert(idx1[i]>=0 && idx1[i]<admat.size());
					assert(idx2[i]>=0 && idx2[i]<admat.size());
					admat[idx1[i]][idx2[i]] = val[i];
				}
				delete dataset_idx1;
				delete dataset_idx2;
				delete dataset_val;
			}
			catch( FileIException not_found_error )
			{
				fprintf(stderr, "Admat not found\n");
			}
		}

		try
		{
			char d_name[1000];
			sprintf(d_name, "%s/transcript_names", group_name);
			DataSet dataset = file->openDataSet(d_name);

    		DataType dtype = dataset.getDataType();
    		assert(H5Tequal(H5Tget_native_type(dtype.getId(), H5T_DIR_DEFAULT), tid1.getId()));

   			DataSpace sid1 = dataset.getSpace();
			int rank = sid1.getSimpleExtentNdims();
			assert(rank==1);
			hsize_t		dims1[rank];
			sid1.getSimpleExtentDims( dims1, NULL);
			
			char* tr_names[dims1[0]];

			dataset.read((void*) tr_names, dtype);

			for(uint i=0; i<dims1[0]; i++) 
			{
				transcript_names.push_back(string(tr_names[i]));
			}
    		dataset.close();
		}
		catch( FileIException not_found_error )
		{
			fprintf(stderr, "transcript_names not found\n");
		}

		try
		{
			// parse transcript paths
			char d_name[1000];
			sprintf(d_name, "%s/transcripts", group_name);
			DataSet* dataset = new DataSet(file->openDataSet(d_name));

			DataType dtype = dataset->getDataType();
			assert(H5Tequal(H5Tget_native_type(dtype.getId(), H5T_DIR_DEFAULT), PredType::NATIVE_INT.getId()));

			DataSpace dataspace = dataset->getSpace();
			assert(dataspace.getSimpleExtentNdims()==2);

			hsize_t dims[2];
			dataspace.getSimpleExtentDims(dims, NULL);
			int num_trans = dims[0];
			int num_seg = dims[1];
			assert(num_seg==segments.size());

			int trans[num_trans][num_seg];
			dataset->read(trans, PredType::NATIVE_INT);

			for (int i=0; i<num_trans; i++)
			{
				vector<int> tmp;
				for (uint j=0; j<num_seg; j++)
				{
					if (trans[i][j]>0)
						tmp.push_back(j);
				}
				transcript_paths.push_back(tmp);
			}
			delete dataset;
		}
		catch( FileIException not_found_error )
		{
			fprintf(stderr, "transcript_paths not found\n");
		}

		delete group;
		delete metadata;
		delete file;
		
	}  // end of try block
	
	// catch failure caused by the H5File operations
	catch( FileIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSet operations
	catch( DataSetIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error )
	{
		error.printError();
		return -1;
	}
	
	return 0;
}
int Bam_Region::write_HDF5(char* filename)
{
    /* First structure  and dataset*/
	typedef struct meta_info {
	int    num_graphs;
	} meta_info;
	
	try
	{
		// suppress printing of exceptions to allow for correct handling
		Exception::dontPrint();

		const H5std_string FILE_NAME(filename);
		H5File* file;
		try 
		{
			// this failes if file does not exist
			file = new H5File(FILE_NAME, H5F_ACC_RDWR);// read and write access to existing file
		}
		catch (Exception error)
		{
			file = new H5File(FILE_NAME, H5F_ACC_TRUNC);// overwrite existing data if file exists
		}

		// read the number of graphs currently stored in the file from the file's 
		// meta data 
		CompType minfo(sizeof(meta_info));
		const H5std_string member_name( "num_graphs" );
		minfo.insertMember( member_name, HOFFSET(meta_info, num_graphs), PredType::NATIVE_INT);
		
		meta_info meta[1];
		DataSet*	metadata;

		const	H5std_string DATASET_NAME( "Graph_meta_info" );
		
		try 
		{  // to determine if the dataset exists in the file
			metadata = new DataSet(file->openDataSet(DATASET_NAME));
			metadata->read( meta, minfo );
		}
		catch( FileIException not_found_error )
		{
			printf("Dataset not found. creating it\n");
			int rank = 1;
			hsize_t dim[] = {1};   /* Dataspace dimensions */
			DataSpace space( rank, dim );
			metadata = new DataSet(file->createDataSet(DATASET_NAME, minfo, space));
			meta[0].num_graphs = 0;
		}
		
		printf("num graphs in file:%i \n", meta[0].num_graphs);

		/*
		* Initialize the data
		*/
		meta[0].num_graphs += 1;
		char group_name[1000];
		sprintf(group_name, "/Graph_%i", meta[0].num_graphs);
		Group* group = new Group(file->createGroup(group_name));

		// Create string datatype
    	StrType tid1(0, H5T_VARIABLE);
    	if(H5T_STRING!=H5Tget_class(tid1.getId()) || !H5Tis_variable_str(tid1.getId()))
       		printf("this is not a variable length string type!!!");
	
		//Write number of graphs to file
		metadata->write( &meta, minfo );

		{
			// write chromosome name to group
			hsize_t		dims1[] = {2};
			int rank = 1;
   			DataSpace sid1(rank, dims1);

			/* Create a dataset */
			char d_name[1000];
			sprintf(d_name, "%s/region_str", group_name);
			DataSet dataset = file->createDataSet(d_name, tid1, sid1);

			char* reg_str[2]; 
			reg_str[0] = get_region_str();
			reg_str[1] = new char[10];
			sprintf(reg_str[1], "%c", strand);
    		/* Write dataset to disk */
    		dataset.write((void*) reg_str, tid1);

			delete[] reg_str[0];
			delete[] reg_str[1];

    		/* Close Dataset */
    		dataset.close();
		}
		if (admat.size()>0)
		{
			// write admat to group
			// first make it sparse
			int len = admat.size();
			vector<int> idx1;
			vector<int> idx2;
			vector<float> val;
			for (int i=0; i<len; i++)
			{
				for (int j=i+1; j<len; j++)
				{
					if (admat[i][j]>=NEIGHBOR)
					{
						idx1.push_back(i);
						idx2.push_back(j);
						val.push_back(admat[i][j]);
					}
				}
			}

			hsize_t dims[1];
			char d_name[1000];
			dims[0] = idx1.size();
			sprintf(d_name, "%s/admat_idx1", group_name);
			DataSpace* s_admat = new DataSpace(1, dims);
			DataSet* d_admat = new DataSet(file->createDataSet(d_name, PredType::NATIVE_INT, *s_admat));
			d_admat->write(&idx1[0], PredType::NATIVE_INT);
			delete s_admat;
			delete d_admat;

			sprintf(d_name, "%s/admat_idx2", group_name);
			dims[0] = idx2.size();
			s_admat = new DataSpace(1, dims);
			d_admat = new DataSet(file->createDataSet(d_name, PredType::NATIVE_INT, *s_admat));
			d_admat->write(&idx2[0], PredType::NATIVE_INT);
			delete s_admat;
			delete d_admat;

			sprintf(d_name, "%s/admat_val", group_name);
			dims[0] = val.size();
			s_admat = new DataSpace(1, dims);
			d_admat = new DataSet(file->createDataSet(d_name, PredType::NATIVE_FLOAT, *s_admat));
			d_admat->write(&val[0], PredType::NATIVE_FLOAT);
			delete s_admat;
			delete d_admat;
		}
		if (segments.size()>0)
		{
			int seg[segments.size()][3];
			for (int i=0; i<segments.size(); i++)
			{
				seg[i][0] = segments[i].first;
				seg[i][1] = segments[i].second;
				seg[i][2] = segments[i].flag;
			}

			int rank = 2;
			hsize_t dims[2];
			char d_name[1000];
			dims[0] = segments.size();
			dims[1] = 3;
			sprintf(d_name, "%s/segments", group_name);
			DataSpace* s_admat = new DataSpace(rank, dims);
			DataSet* d_admat = new DataSet(file->createDataSet(d_name, PredType::NATIVE_INT, *s_admat));
			d_admat->write(seg, PredType::NATIVE_INT);
			delete s_admat;
			delete d_admat;

		}
		if (transcript_names.size()>0)
		{
			int num = transcript_names.size();
	 		char* wdata[num];
			for (int i=0; i<num; i++)
			{
				wdata[i] = new char[transcript_names[i].length()+1];
				sprintf(wdata[i], "%s", transcript_names[i].c_str());
			}

    		hsize_t		dims1[] = {num};
			int rank = 1;
   			DataSpace sid1(rank, dims1);

			/* Create a dataset */
			char d_name[1000];
			sprintf(d_name, "%s/transcript_names", group_name);
			DataSet dataset = file->createDataSet(d_name, tid1, sid1);

    		/* Write dataset to disk */
    		dataset.write((void*) wdata, tid1);

			for (int i=0; i<num; i++)
				delete[] wdata[i];
    		/* Close Dataset */
    		dataset.close();
	
		}
		if (transcript_paths.size()>0)
		{
			int num_seg = segments.size();
			int num_trans = transcript_paths.size();
			int trans[num_trans][num_seg];

			memset(trans, 0, num_trans*num_seg*sizeof(int));
			for (int i=0; i<num_trans; i++)
			{
				for (uint j=0; j<transcript_paths[i].size(); j++)
				{
					int seg = transcript_paths[i][j];
					assert(seg<num_seg);
					trans[i][seg]=1;
				}
			}
			int rank = 2;
			hsize_t dims[2];
			char d_name[1000];
			dims[0] = num_trans;
			dims[1] = num_seg;
			sprintf(d_name, "%s/transcripts", group_name);
			DataSpace* s_admat = new DataSpace(rank, dims);
			DataSet* d_admat = new DataSet(file->createDataSet(d_name, PredType::NATIVE_INT, *s_admat));
			d_admat->write(trans, PredType::NATIVE_INT);
			delete s_admat;
			delete d_admat;
		}
		
		/*
		* Release resources
		*/
		delete group;
		delete metadata;
		delete file;
		
	}  // end of try block
	
	// catch failure caused by the H5File operations
	catch( FileIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSet operations
	catch( DataSetIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error )
	{
		error.printError();
		return -1;
	}
	
	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error )
	{
		error.printError();
		return -1;
	}
	
	return 0;
}
int Bam_Region::write_binary(std::ofstream* ofs)
{
	// region meta info
	ofs->write((char *) &start, sizeof(int));
	ofs->write((char *) &stop, sizeof(int));
	int len = strlen(chr);
	assert(len>0);
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) chr, len);
	ofs->write((char *) &strand, sizeof(char));

	// splice graph
	len = segments.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &segments[0], len*sizeof(segment));
	if (seg_cov.size()==0)
		compute_seg_cov();
	//if (len==1 && seg_cov[0]<1e-3)
	//{
	//	printf("seg_cov for single segment locus: %.2f\n", seg_cov[0]);
	//	printf("number of reads (filtered): %i\n", (int)reads.size());
	//}
	ofs->write((char *) &seg_cov[0], len*sizeof(float));
	
	//splice graph admat
	len = admat.size();
	ofs->write((char *) &len, sizeof(int));
	vector<int> sp;
	vector<float> val;
	for (int i=0; i<len; i++)
		for (int j=i+1; j<len; j++)
			if (admat[i][j]>NEIGHBOR)
			{
				sp.push_back(i);
				sp.push_back(j);
				val.push_back(admat[i][j]);
			}
	len = sp.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &sp[0], len*sizeof(int));
	len = val.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &val[0], len*sizeof(float));


	// transcript paths
	int num_trans = transcript_paths.size();
	ofs->write((char *) &num_trans, sizeof(int));
	for (int i=0; i<num_trans; i++)
	{
		len = transcript_paths[i].size();
		ofs->write((char *) &len, sizeof(int));
		if (len>0)
		{
			ofs->write((char *) &transcript_paths[i][0], len*sizeof(int));
		}
	}

	// transcript names
	assert(transcript_paths.size()==transcript_names.size());
	for (int i=0; i<num_trans; i++)
	{
		int len = transcript_names[i].length();
		ofs->write((char *) &len, sizeof(int));
		ofs->write(transcript_names[i].c_str(), len*sizeof(char));
	}

	// pair matrix
	len = pair_mat.size();
	ofs->write((char *) &len, sizeof(int));
	if (len==0)
		return -1;
	sp.clear();
	vector<int> pair_cnt;
	for (int i=0; i<len; i++)
		for (int j=i+1; j<len; j++)
			if (pair_mat[i][j]>0)
			{
				sp.push_back(i);
				sp.push_back(j);
				pair_cnt.push_back(pair_mat[i][j]);
			}
	len = sp.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &sp[0], len*sizeof(int));
	len = pair_cnt.size();
	ofs->write((char *) &len, sizeof(int));
	ofs->write((char *) &pair_cnt[0], len*sizeof(int));

	return 1;
}
int Bam_Region::read_binary(std::ifstream* ifs)
{
	// region meta info
	ifs->read((char *) &start, sizeof(int));
	ifs->read((char *) &stop, sizeof(int));
	int len = 0;
	ifs->read((char *) &len, sizeof(int));
	if (len<=0)
	{
		fprintf(stderr, "read_binary: could not read from file: len:%i\n", len);
		fprintf(stderr, "read_binary: start:%i stop:%i\n", start, stop);
		return -1;
	}
	chr = new char[len+1];
	ifs->read((char *) chr, len);
	chr[len] = '\0';
	ifs->read((char *) &strand, sizeof(char));

	//printf("%s:%i..%i %c\n", chr, start, stop, strand);
	// splice graph
	ifs->read((char *) &len, sizeof(int));
	assert(len>0);
	segment seg[len];
	ifs->read((char *) seg, len*sizeof(segment));
	for (int i=0; i<len; i++)
	{
		segments.push_back(seg[i]);
	}
	seg_cov.resize(len);
	ifs->read((char *) &seg_cov[0], len*sizeof(float));

	//splice graph admat
	ifs->read((char *) &len, sizeof(int));
	assert(len>0);
	assert(admat.size()==0);
	for (int i=0; i<len; i++)
	{
		vector<float> row(len, -2);
		admat.push_back(row);
	}
	ifs->read((char *) &len, sizeof(int));
	assert(len>0);
	int sp[len];
	ifs->read((char *) sp, len*sizeof(int));
	int len2 = 0;
	ifs->read((char *) &len2, sizeof(int));
	assert(len2*2==len);
	float val[len2];
	ifs->read((char *) val, len2*sizeof(float));
	for (int i=0; i<len2; i++)
	{
		int j = sp[2*i];
		int k = sp[2*i+1];
		admat[j][k] = val[i];
		admat[k][j] = val[i];
	}
	// restore connections of neighboring segments
	for (uint i=1; i<segments.size(); i++)
	{
		if (segments[i-1].second+1==segments[i].first)
		{
			admat[i][i+1] = NEIGHBOR;
		}
	}
		
	fd_out = stdout;

	// transcript paths
	int num_trans = 0;
	ifs->read((char *) &num_trans, sizeof(int));
	for (int i=0; i<num_trans; i++)
	{
		ifs->read((char *) &len, sizeof(int));
		if (len>0)
		{
			vector<int> path(len, 0);
			ifs->read((char *) &path[0], len*sizeof(int));

			for (uint j=0; j<path.size(); j++)
			{
				assert(path[j]<(int)segments.size());
				if (j>0)
					assert(admat[path[j-1]+1][path[j]+1]>=NEIGHBOR);
			}
			transcript_paths.push_back(path);
		}
	}

	for (int i=0; i<num_trans; i++)
	{
		int len;
		ifs->read((char *) &len, sizeof(int));
		char name[len+1];
		ifs->read(name, len*sizeof(char));
		string sname(name, len);
		transcript_names.push_back(sname);
	}

	//pair_mat
	ifs->read((char *) &len, sizeof(int));
	if (len==0)
		return 0;

	pair_mat.clear();
	for (uint i=0; i<segments.size(); i++)
	{
		vector<int> row(segments.size(), 0);
		pair_mat.push_back(row);
	}
	int num_pairs = 0;
	ifs->read((char *) &num_pairs, sizeof(int));
	int pairs[num_pairs];
	ifs->read((char *) pairs, num_pairs*sizeof(int));
	ifs->read((char *) &len, sizeof(int));
	int* pair_cnt = new int[len];
	ifs->read((char *) pair_cnt, len*sizeof(int));

	assert(num_pairs==2*len);
	//for (int i=0; i<len; i+=2)
	//{
	//	int j = pairs[i]; 
	//	int k = pairs[i+1];
	//	if (j<0 || j >= pair_mat.size() || k<0 || k>=pair_mat.size())
	//	{
	//		fprintf(stderr, "Error: j: %i k: %i, pair_mat.size(): %i\n", j, k, (int) pair_mat.size());
	//		break;
	//	}
	//	pair_mat[j][k] = pair_cnt[i/2];
	//}
	delete[] pair_cnt;
	return 0;
}
void Bam_Region::compute_intron_list()
{
	intron_list.clear();
	unique_introns.clear();
   	intron_counts.clear();

	vector<int> introns;
	for (uint i=0; i<reads.size(); i++)
	{
		reads[i]->get_introns(&introns);
    }
	for (uint i=0; i<introns.size(); i+=2)
	{
		segment intr(introns[i], introns[i+1]);
		intron_list.push_back(intr);
	}
	fprintf(fd_out, "found %lu introns\n", intron_list.size());
	
	// sort by intron start
	sort(intron_list.begin(), intron_list.end(), compare_second);

	if (intron_list.size()>0)
	{
    	unique_introns.push_back(intron_list[0]);
		intron_counts.push_back(1);
    	for (uint i=1; i<intron_list.size(); i++)
    	{
			//printf("intron: %i->%i\n", intron_list[i].first, intron_list[i].second);
    	    if (intron_list[i].first!=intron_list[i-1].first || intron_list[i].second!=intron_list[i-1].second)
    	    {
				//printf("uintron: %i->%i\n", unique_introns.back().first, unique_introns.back().second);
    	       	unique_introns.push_back(intron_list[i]);
				intron_counts.push_back(1);
    	    }
			else
			{
    	    	intron_counts.back()++;
			}
    	}
	}

	//for (int i=1; i<unique_introns.size(); i++)
	//	assert(unique_introns[i].first!=unique_introns[i-1].first || unique_introns[i].second!=unique_introns[i-1].second);

	//vector<segment*> same_start;
	//for (int i=0; i<intron_list.size(); i++)
	//{
	//	if (i<intron_list.size()-1 && intron_list[i].first==intron_list[i+1].first)
	//	{
	//		// collect introns with same start
	//		same_start.push_back(&intron_list[i]);
	//	}
	//	else
	//	{
	//		same_start.push_back(&intron_list[i]);
	//		sort(same_start.begin(), same_start.end(), compare_second);
	//		unique_introns.push_back(*(same_start[0]));
	//		intron_counts.push_back(1);
	//		for (int j=1; j<same_start.size(); j++)
	//		{
	//			if (same_start[j]->second!=same_start[j-1]->second)
	//			{
	//				unique_introns.push_back(*(same_start[j]));
	//				intron_counts.push_back(1);
	//			}
	//			else
	//			{
	//				//printf("intron_counts before: %i ", (*(intron_counts.back())));
	//				(intron_counts.back())++; 
	//				//printf("after: %i \n", (*intron_counts.back()));
	//			}
	//		}
	//		same_start.clear();
	//	}
	//}
	fprintf(fd_out, "found %lu unique introns\n", unique_introns.size());

#ifdef READ_DEBUG
	//check unique list
	printf("DEBUG: Check 'unique introns'-list contains all introns\n");
	uint idx = 0;
	for (int i=0; i<intron_list.size(); i++)
	{
		bool match = false;
		while (unique_introns[idx].first<intron_list[i].first && idx<unique_introns.size())
			idx++;
		uint idx2 = idx;
		while (unique_introns[idx2].first==intron_list[i].first && idx2<unique_introns.size())
		{
			if (unique_introns[idx2].second==intron_list[i].second)
				match = true;
			idx2++;
		}
		assert(match);
	}
#endif

	int sum = 0;
	for (uint i=0; i<intron_counts.size(); i++)
		sum+=intron_counts[i]; 
	//if (sum!=intron_list.size())
	//{
	//	printf("Error %i!=%i\n", sum, (int)intron_list.size());
	//	for (int i=0; i<intron_list.size(); i++)
	//	{
	//		printf("intron: %i->%i\n", intron_list[i].first, intron_list[i].second);
	//	}
	//	for (int i=0; i<unique_introns.size(); i++)
	//	{
	//		printf("uintron: %i->%i\t%i\n", unique_introns[i].first, unique_introns[i].second, intron_counts[i]);
	//	}
	//}
	assert(sum==(int)intron_list.size());
}

