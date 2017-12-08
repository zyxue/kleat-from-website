#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <cstdlib> //for atoi()

using namespace std;

namespace opt {
	static bool sense;
	static unsigned support;
	static unsigned pet_support;
}

const static int pad = 100;

struct overlapping_range
{
	pair<int, int> range;

	overlapping_range(int start, int end) {
		range = make_pair(start, end);
	}

	overlapping_range(pair<int, int> range) : range(range) {}

	bool isOverlaping(overlapping_range o_range) {
		return range.first <= o_range.range.second + pad &&
			range.second >= o_range.range.first - pad;
	}

	bool isOverlaping(int val) {
		return range.first <= val + pad &&
			range.second >= val - pad;
	}

	overlapping_range overlap(overlapping_range o_range) {
		pair<int, int> new_range = range;
		if (isOverlaping(o_range))
			new_range = make_pair(
					min(range.first, o_range.range.first),
					max(range.second, o_range.range.second));
		return overlapping_range(new_range);
	}

	bool operator==(const overlapping_range o_range) const {
		return range.first == o_range.range.first &&
			range.second == o_range.range.second;
	}

	bool operator<(const overlapping_range o_range) const {
		return range.first < o_range.range.first;
	}
};

struct BedGraph {
	typedef vector<int> Covs;
	typedef string Chrom;
	map<Chrom, Covs> values;
	map<Chrom, vector<bool> > peaks;
	map<Chrom, vector<pair<overlapping_range, bool> > > overlaps;

	int max_element(Chrom& chr, int start, int end)
	{
		Covs::iterator it = values[chr].begin();
		int size = values[chr].size();
		if (start > size)
			return 0;
		if (end > size)
			end = size;
		return *std::max_element(it + start, it + end);
	}

	bool isOverlaping(Chrom chr, int start) {
		for (vector<pair<overlapping_range, bool> >::iterator it = overlaps[chr].begin();
				it < overlaps[chr].end(); it++) {
			if (it->first.isOverlaping(start)) {
				it->second = true;
				return true;
			}
		}
		return false;
	}

	void setPeaks(int thresh) {
		cerr << "Finding peaks...\n";
		for (map<Chrom, Covs>::iterator mit = values.begin();
				mit != values.end(); mit++) {
			bool prev = 0;
			pair<int, int> curr_range;
			Chrom chr = mit->first;
			Covs c = mit->second;
//			vector<bool> v(c.size());
			vector<overlapping_range> ranges;
			ranges.reserve(c.size());
			for (int i = 0; i < c.size(); i++) {
				bool yup = c[i] >= thresh;
//				if (yup)
//					v[i] = 1;
				if (prev ^ yup) {
					if (!prev)
						curr_range.first = i;
					else {
						curr_range.second = i;
						ranges.push_back(overlapping_range(curr_range));
					}
				}
				prev = yup;
			}
			sort(ranges.begin(), ranges.end());
			vector<pair<overlapping_range, bool> > nranges;
			for (unsigned i = 1; i < ranges.size(); i++) {
				overlapping_range r = ranges[i-1];
				while (r.isOverlaping(ranges[i]) && i < ranges.size()) {
					r = r.overlap(ranges[i]);
					i++;
				}
				nranges.push_back(make_pair(r, false));
			}
//			peaks.insert(make_pair(chr, v));
			overlaps.insert(make_pair(chr, nranges));
			cerr << "Found peaks for " << chr << ".\n";
			cerr << "sizeof nranges: " << nranges.size()
				<< " sizeof ranges: " << ranges.size()
				<< " sizeof c: " << c.size() << '\n';
		}
		cerr << "Found the peaks.\n";
	}

	friend istream& operator>>(istream& in, BedGraph& bg)
	{
		int i = 0;
		string s;
		while (getline(in, s)) {
			stringstream ss(s);
			if (ss.peek() == '#')
				continue;

			string chr;
			int start, end, val;
			ss >> chr >> start >> end >> val;
			if (end > bg.values[chr].size())
				bg.values[chr].resize(end);

			for (int x = start; x < end; x++)
				bg.values[chr][x] += val;

			if (++i % 100000 == 0)
				cerr << "loaded " << i << " entries.\n";
		}
		return in;
	}
	
	friend ostream& operator<<(ostream& out, BedGraph& bg) {
		for (map<string, vector<pair<overlapping_range, bool> > >::iterator it = bg.overlaps.begin();
				it != bg.overlaps.end(); it++) {
			string chr = it->first;
			for (vector<pair<overlapping_range, bool> >::iterator vit = it->second.begin();
					vit != it->second.end(); vit++) {
				//int d = distance(it->second.begin(), vit);
				//if (!bg.peaks[chr][d])
				//	continue;
				out << chr << '\t'
					<< vit->first.range.first << '\t'
					<< vit->first.range.second << '\t'
					<< "TP=" << vit->second << '\n';
				//out << chr << '\t'
				//	<< d <<'\t'
				//	<< bg.values[chr][d] << '\t'
				//	<< *vit << '\n';
			}
		}
		return out;
	}
};



void handleEvents(const char* eventsPath, BedGraph& bg)
{
	string line;
	ifstream in(eventsPath);
	while(getline(in, line)) {
		stringstream ss(line);
		string chr;
		ss >> chr;
		if (chr == "track")
			continue;
		int start, end;
		float support;
		char strand;
		string gene;
		ss >> start >> end >> support >> strand >> gene;
		if (support < opt::support)
			continue;
		if ((strand == '+') ^ opt::sense)
			continue;
		bool o = bg.isOverlaping(chr, start);

		if (opt::sense)
			end += 100;
		else
			start -= 100;
		start = start < 0 ? 0 : start;
		int max_e = bg.max_element(chr, start, end);
		cout << chr << '\t'
			<< start << '\t'
			<< end << '\t'
			<< max_e << '\t'
			<< support << '\t'
			<< o << '\t'
			<< gene << '\n';
	}
	assert(in.eof());
}

int main(int argc, char** argv)
{
	if (argc != 6) {
		cerr << "Usage: prog <RNA-PET.bedgraph> <polya.bg> [+-] <support> <pet_support>\n";
		exit(EXIT_FAILURE);
	}
	ifstream bgfile(argv[argc-5]);
	BedGraph bg;
	bgfile >> bg;

	opt::sense = string(argv[argc-3]) == "+";
	assert(opt::sense || string(argv[argc-3]) == "-");
	opt::support = atoi(argv[argc-2]);
	opt::pet_support = atoi(argv[argc-1]);
	bg.setPeaks(opt::pet_support);

	//cout << bg;

	handleEvents(argv[argc-4], bg);

	cout << bg;

	return 0;
}
