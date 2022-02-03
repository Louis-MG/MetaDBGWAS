/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/** \file Histogram.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Histogram feature
 */

#ifndef _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_
#define _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/collections/api/Bag.hpp>
#include <gatb/tools/misc/api/IHistogram.hpp>
#include <string>
#include <iostream>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace tools     {
namespace misc      {
namespace impl      {
/********************************************************************************/

/** \brief Basic implementation of the IHistogram interface.
 *
 * This implementation is the one actually used by SortingCountAlgorithm.
 */
class Histogram : public IHistogram, public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] length : maximum value for the X axis
     * \param[in] bag : bag where the values can be saved. */
    Histogram (size_t length)
        : _length(length), _cutoff(0), _nbsolids(0), _firstPeak(0),
          _histogram(0), _histogram_smoothed(0)
    {
        _histogram = (Entry*) CALLOC (_length + 1, sizeof (Entry));
        memset (_histogram, 0, sizeof(Entry)*(_length + 1));

		_histogram_smoothed = (Entry*) CALLOC (_length + 1, sizeof (Entry));
        memset (_histogram_smoothed, 0, sizeof(Entry)*(_length + 1));
		
        for (size_t i=0; i<_length+1; i++)
        {
            _histogram[i].index     = i;
			_histogram_smoothed[i].index     = i;

            _histogram[i].abundance = 0;
        }
    }

    /** Destructor */
    virtual ~Histogram ()
    {
        FREE (_histogram);
        FREE (_histogram_smoothed);
    }

    /** \copydoc IHistogram::inc */
    void inc (u_int16_t index)  { _histogram [(index >= _length) ? _length : index].abundance ++; }

    /** \copydoc IHistogram::save */
    void save (tools::storage::impl::Group& group);

    /** \copydoc IHistogram::compute_threshold */
	void compute_threshold (int min_auto_threshold) ;  //min_auto_threshold = prevents the auto_cutoff from being below this value. Default =3
	
    /** \copydoc IHistogram::get_solid_cutoff */
	u_int16_t get_solid_cutoff ()  { return _cutoff; }

    /** \copydoc IHistogram::get_nbsolids_auto */
	u_int64_t get_nbsolids_auto ()  { return _nbsolids; }

    /** \copydoc IHistogram::get_first_peak */
	u_int16_t get_first_peak ()  { return _firstPeak; }

    /** \copydoc IHistogram::getLength */
    size_t getLength() { return _length; }

    /** \copydoc IHistogram::get */
    u_int64_t& get (u_int16_t idx)  { return _histogram[idx].abundance; }

private:

    size_t    _length;
	u_int16_t _cutoff;
	u_int64_t _nbsolids;
    u_int16_t _firstPeak;
	
    Entry*  _histogram;
	Entry*  _histogram_smoothed;
};

/********************************************************************************/

/** \brief Null implementation of the IHistogram interface.
 */
class HistogramNull : public IHistogram, public system::SmartPointer
{
public:

    /** \copydoc IHistogram::inc */
    void inc (u_int16_t index) {}

    /** \copydoc IHistogram::save */
    void save (tools::storage::impl::Group& group)  {}
	
    /** \copydoc IHistogram::get_solid_cutoff */
	u_int16_t get_solid_cutoff  () { return 0; }

    /** \copydoc IHistogram::get_nbsolids_auto */
	u_int64_t get_nbsolids_auto () { return 0; }

    /** \copydoc IHistogram::get_first_peak */
	u_int16_t get_first_peak    () { return 0; }

    /** \copydoc IHistogram::compute_threshold */
	void compute_threshold (int min_auto_threshold) { }

    /** \copydoc IHistogram::getLength */
    size_t getLength() { return 0; }

    /** \copydoc IHistogram::get */
    u_int64_t& get (u_int16_t idx)  { static u_int64_t foo; return foo; }
};

/********************************************************************************/

/** \brief Cached implementation of the IHistogram interface.
 *
 * This implementation is a Proxy design pattern. It allows to modify a IHistogram instance
 * by several threads at the same time. Actually, each thread has a local copy and at
 * the end, all the local copies are merged into the referred instance.
 * */
class HistogramCache : public IHistogram, public system::SmartPointer
{
public:

    /** Constructor.
     * \param[in] ref : the referred instance.
     * \param[in] synchro : used for synchronization */
    HistogramCache (IHistogram* ref, system::ISynchronizer* synchro=0)
        : _ref(0), _synchro(synchro), _localHisto(ref ? ref->getLength() : 0) {  setRef(ref); }

    /** Destructor. */
    ~HistogramCache()
    {
        system::LocalSynchronizer ls (_synchro);
        for (size_t cc=1; cc<_localHisto.getLength(); cc++)  {  _ref->get(cc) += _localHisto.get(cc);  }
        setRef (0);
    }

    /** \copydoc IHistogram::inc */
    void inc (u_int16_t index)  { _localHisto.inc (index); }

    /** \copydoc IHistogram::save */
    void save (tools::storage::impl::Group& group)  { return _ref->save(group); }

    /** \copydoc IHistogram::compute_threshold */
	void compute_threshold (int min_auto_threshold) { return _ref->compute_threshold(min_auto_threshold); }

    /** \copydoc IHistogram::get_solid_cutoff */
	u_int16_t get_solid_cutoff () {return _ref->get_solid_cutoff();}
	
    /** \copydoc IHistogram::get_nbsolids_auto */
	u_int64_t get_nbsolids_auto () {return _ref->get_nbsolids_auto();}

    /** \copydoc IHistogram::get_first_peak */
    u_int16_t get_first_peak () { return _ref->get_first_peak(); }

    /** \copydoc IHistogram::getLength */
    size_t getLength() { return _localHisto.getLength(); }

    /** \copydoc IHistogram::get */
    u_int64_t& get (u_int16_t idx)  { return _localHisto.get(idx); }

private:

    IHistogram* _ref;
    void setRef (IHistogram* ref)  { SP_SETATTR(ref); }

    system::ISynchronizer* _synchro;
    Histogram              _localHisto;
};

/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_HISTOGRAM_HPP_ */
