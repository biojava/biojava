package org.biojava3.genome.parsers.gff;


/**
 * A location on a sequence.
 * A location is a contiguous range of indices, with a single start and end point. 
 * <br><br>
 * Internally, location indices are stored in Java "half-open" format: the start is the (origin 0) index of
 * the first symbol in the range; the end is the origin 0 index of the first symbol PAST the
 * end of the range, so that end - start == length.
 * <br><br>
 * Location objects, once constructed, cannot be changed. Instead, all methods return a new
 * location. This allows the use of "method chaining" to implement a particular calculation.
 * For example, consider the chained statement "loc.prefix( 100 ).suffix( 10 )", 
 * which first applies the prefix method to
 * the variable named loc, and then the suffix method to the result. 
 * Together, the chained operations create a new location object of length 10
 * whose start is the index of the 90th symbol.
 * Here's another example. This one returns a location object holding the coordinates of the intron between 
 * the first exon (location exon1) and 
 * the second exon (location exon2) on a sequence (seq): "seq.prefix( exon2 ).suffix( exon1 )"
 * <br><br>
 * About the negative (reverse) strand: The location object stores reverse strand locations as
 * negative indices. For example, the positive strand location from index 12 to index 97 is
 * on the opposite side as index -97 (start) to index -12 (end). Note that the larger index is 
 * always downstream from the smaller index, (i.e. start &lt;= end, regardless of strand). 
 * Obviously this representation makes it trivial
 * to convert a location from one strand to the other.
 * <br><br>
 * Additional points regarding the use of locations on opposite strands:
 *<br>
 * (1) Opposite strand locations cannot be compared, eg isBefore() will
 * throw an exception.<br>
 * (2) Containment queries ( eg overlaps(), contains() ) also throw exceptions.
 *<br>
 * (3) The plus() method will map a location to its positive strand equivalent; use it on both args 
 * before calling, for example the intersection() method,
 * if your code needs to be indifferent to strand.
 *<br><br>
 * Exceptions and how they are (typically) used:
 *<br>
 * IllegalArgumentException - the location given as a parameter is not on the same strand as the location.
 *<br>
 * IndexOutOfBoundsException - often means the operation caused the location to span the origin, ie
 * be partially on positive and partially on negative strand.
 *<br>
 * @author Hanno Hinsch
 */
public class Location implements Iterable<Location>
{
	
	private int mStart;
	private int mEnd;

//	/**
//	 */
//	private Location()
//	{
//	}
	
	/**
	 * Construct new location from coordinates.
	 * See package description of coordinate format.
	 * @param start Origin 0 index of first symbol.
	 * @param end Origin 0 index of last symbol + 1.
	 * @throws IllegalArgumentException End is not after start, or location spans the origin
	 */
	public Location( int start, int end )
	{
		mStart= start;
		mEnd= end;
		
		if( !isHealthy() )
		{
			throw new IllegalArgumentException( "Improper location parameters: (" + start + "," + end + ")" );
		}
		
	}
	
	/**
	 * Clone other location.
	 * @param other The location to clone.
	 */	
	public Location( Location other )
	{
		mStart= other.mStart;
		mEnd= other.mEnd;
		
		assert isHealthy(): toString();
	}

        public int getBegin(){
            if(isNegative())
                return mEnd;
            else
                return mStart;
        }

        public int getEnd(){
            if(isNegative())
                return mStart;
            else
                return mEnd;
        }


	/**
	 * Create location from "biocoordinates", as in GFF file. In biocoordinates,
	 * the start index of a range is represented in origin 1 (ie the very first index is 1, not 0),
	 * and end= start + length - 1. 
	 *
	 * @param start Origin 1 index of first symbol.
	 * @param end Origin 1 index of last symbol.
	 * @param strand '+' or '-' or '.' ('.' is interpreted as '+').
	 * @return Created location.
	 * @throws IllegalArgumentException strand must be '+', '-' or '.'
	 */
	public static Location fromBio( int start, int end, char strand )
	{
		int s= start - 1;
		int e= end;
		
		if( !( strand == '-' || strand == '+' || strand == '.' ))
		{
			throw new IllegalArgumentException( "Strand must be '+', '-', or '.'" );
		}
		
		if( strand == '-' )
		{
			//negate
			s= - end;
			e= - ( start - 1);
		}
		
		return new Location( s, e );
	}
	
	/**
	 * Create a location from MAF file coordinates, which represent negative
	 * strand locations as the distance from the end of the sequence.
	 *
	 * @param start Origin 1 index of first symbol.
	 * @param length Number of symbols in range.
	 * @param strand '+' or '-' or '.' ('.' is interpreted as '+').
	 * @param totalLength Total number of symbols in sequence.
	 * @throws IllegalArgumentException Strand must be '+', '-', '.'
	 * 
	 */
	 public static Location fromBioExt( int start, int length, char strand, int totalLength )
	 {
		int s= start;
		int e= s + length;
		
		if( !( strand == '-' || strand == '+' || strand == '.' ))
		{
			throw new IllegalArgumentException( "Strand must be '+', '-', or '.'" );
		}
		
		if( strand == '-' )
		{
			s= s - totalLength;
			e= e - totalLength;
		}
		
		return new Location( s, e );
	 }
	 
	/**
	 * Get character representation of strand.
	 *
	 * @return '+' or '-'
	 */
	public char bioStrand()
	{
		return ( isNegative() )?'-':'+';
	}
	
	/**
	 * Get start index, in biocoordinates.
	 *
	 * @return The origin 1 index of the first symbol in location.
	 */
	public int bioStart()
	{
		return plus().start() + 1;
	}
	
	/**
	 * Get end index, in biocoordinates.
	 *
	 * @return The origin 1 index of the final symbol in location.
	 */
	public int bioEnd()
	{
		return plus().end();
	}
	


	/**
	 * Return location that is in same position on plus strand. If location is already
	 * on plus strand, just return the location unchanged.
	 *
	 * @return Location on plus strand.
	 */
	public Location plus()
	{
		if( isNegative() )
		{
			return opposite();
		}
		else
		{
			return this;
		}
	}
	
	/**
	 * Return location that is in same position on negative strand. If location is already
	 * on negative strand, just return the location unchanged.
	 *
	 * @return Location on negative strand.
	 */
	public Location minus()
	{
		if( isNegative() )
		{
			return this;
		}
		else
		{
			return opposite();
		}
	}

	
	/**
	*  Return the union. 
	* <br>
	*
	 * @param other The location to join.
	 * @return The union is a range that starts at the lesser of the two starting indices and ends at the
	 * greater of the two ends.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	*/
	public Location union( Location other )
	{
		
		if( !isSameStrand( other ))
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
		else
		{
			int start= (other.mStart < mStart)? other.mStart: mStart;
			int end= (other.mEnd > mEnd)? other.mEnd: mEnd;
			
			return new Location( start, end );
		}

	}

	/**
	* Return the intersection, or null if no overlap.  
	*	 
	 * @param other The location to intersect.
	 * @return The maximal location that is contained by both. Returns null if no overlap!
	 * @throws IllegalArgumentException Locations are on opposite strands.
	*/
	public Location intersection( Location other )
	{
		if( isSameStrand( other ) )
		{
			int start= ( mStart > other.mStart)? mStart: other.mStart;	//pick larger
			int end= ( mEnd < other.mEnd)? mEnd: other.mEnd;	//pick smaller
				
			if( start <= end )
			{
				return new Location( start, end );
			}
			else
			{
				return null;
			}			
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );				
		}
		
	}
	
	/**	
	 * Get starting index (origin 0).
	 *
	 * @return The start index.
	 */
	public int start()
	{
		return mStart;
	}
	
	/**	
	 * Get the ending index.
	 *
	 * @return The index of last symbol + 1 (remember Java half-open coordinates).
	 */
	public int end()
	{
		return mEnd;
	}
	
	/**	
	 * Get length of range.
	 *
	 * @return The length of the range (end - start).
	 */
	public int length() 
	{ 
		return mEnd - mStart; 
	}
	
	
	/**
	 * Enable a "sliding window" iteration over a location 
	 * to use with Java's "for" loop construct.
	 * The returned helper object implements the Iterable interface; the windowSize and increment semantics are implemented
	 * by an underlying LocIterator.
	 * <br><br>
	 * For example, given a location variable "loc":
	 *<br>
<pre>
	//use window size of 3 and increment of +3
	for( Location temp: loc.window( 3, 3 ))
	{
	//at each iteration, temp will be the location of the next 3 symbols
	}
</pre>
	 *
	 * @param windowSize The number of symbols to get on each iteration.
	 * @param increment The direction and number of symbols to advance at each iteration.
	 * @return An anonymous iterable object to use with Java's for( ... ) loop construct.
	 */
	public Iterable<Location> window( final int windowSize, final int increment )
	{
		final Location loc= this;
		
		//return iterable anonymous inner class
		return new Iterable<Location> ()
			{
				public LocIterator iterator()
				{
					return new LocIterator( loc, windowSize, increment );
				}
				
			};
	}
	
	/**	
	 * Create a location iterator over this location with a window size of 1 and
	 * an increment of +1 (successive symbols from start to end).
	 *
	 * @return An iterator over a Location (a LocIterator object).
	 */
	public LocIterator iterator()
	{
		return new LocIterator( this, 1, 1 );
	}
	
	/**	
	 * Create a location iterator over this location, 
	 * using specified window size and increment.
	 *
	 * @param windowSize The number of symbols to get on each iteration.
	 * @param increment The direction and number of symbols to advance at each iteration.
	 * @return An iterator over a Location (a LocIterator object).
	 */
	public LocIterator iterator( int windowSize, int increment )
	{
		return new LocIterator( this, windowSize, increment );
	}
	
			
	/**
	 * The part of this location before the specified position. If position is negative,
	 * count backwards from the end.
	 * <br><br>
	 * For position >= 0, return Location( start, start + position ).
	 * <br>
	 * For position < 0, return Location( start, end + position ).
	 * <br>
	 * @return New location from start of this location to directly before position.
	 * @param position Where the prefix ends.
	 * @throws IndexOutOfBoundsException Specified prefix is longer than location.
	 */
	public Location prefix( int position )
	{
		int end;
		if( position >= 0 )
		{
			if( (mStart + position <= mEnd) )
			{
				end= mStart + position;
			}
			else
			{
				throw new IndexOutOfBoundsException( "Specified prefix longer than location." );
			}
		}
		else
		{
			if( (mEnd + position > mStart))
			{
				end= mEnd + position;
			}
			else
			{
				throw new IndexOutOfBoundsException( "Specified prefix longer than location." );
			}
		}
		
		return new Location( mStart, end );		
	}
	
	
	/**
	 * The part of this location after the specified position. If position is negative, count backwards
	 * from the end.
	 * <br><br>
	 * For position >= 0, return Location( start + position, end ).
	 * <br>
	 * For position < 0, return Location( end - position, end ).
	 * <br>
	 * @return New location from position to end of this location.
	 * @param position Where the suffix starts.
	 * @throws IndexOutOfBoundsException Specified suffix is longer than location.
	 */
	public Location suffix( int position )
	{
		int start;
		if( position >= 0 )
		{
			if( mStart + position <= mEnd ) // Scooter willis when 60 + 60 = 120 no remainder
			{
				start= mStart + position;
			}
			else
			{
				throw new IndexOutOfBoundsException( "Specified suffix longer than location." );
			}
		}
		else
		{
			if( mEnd + position >= mStart )
			{
				start= mEnd + position;
			}
			else
			{
				throw new IndexOutOfBoundsException( "Specified suffix longer than location." );
			}
		}
		
		return new Location( start, mEnd );
	}
	
	/** 
	 * The part of this location before the other location (not inclusive).
	 *
	 * @param other The other location.
	 * @return The part of this location before the other location.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 * @throws IndexOutOfBoundsException This location does not contain other location.
	 */
	public Location prefix( Location other )
	{
		
		if( isSameStrand( other ) )
		{
			if( other.mStart >= mStart )
			{
				return new Location( mStart, (other.mStart < mEnd)? other.mStart: mEnd );
			}
			else
			{
				//other is out of bounds -- no prefix
				throw new IndexOutOfBoundsException( "Specified location not within this location." );
			}
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**
	 * The part of this location after the other location (not inclusive). 
	 *
	 * @param other The other location.
	 * @return The part of this location after the other location.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 * @throws IndexOutOfBoundsException This location does not contain other location.
	 */
	public Location suffix( Location other )
	{
		if( isSameStrand( other ))
		{
			if( other.mEnd <= mEnd )
			{ 
				return new Location( (other.mEnd > mStart)? other.mEnd: mStart, mEnd );
			}
			else
			{
				//other is out of bounds -- no suffix
				throw new IndexOutOfBoundsException( "Specified location not within this location." );
			}
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
			
	}
	
	/**
	 * Return the adjacent location of specified length directly upstream of this location.
	 *
	 * @return Upstream location.
	 * @param length The length of the upstream location.
	 * @throws IndexOutOfBoundsException Specified length causes crossing of origin.
	 */
	public Location upstream( int length )
	{
		if( length < 0 )
		{
			throw new IllegalArgumentException( "Parameter must be >= 0; is=" + length );
		}
		
		if( Math.signum( mStart - length) == Math.signum( mStart ) || 0 == Math.signum( mStart - length ) )
		{
			return new Location(mStart - length, mStart );
		}
		else
		{
			throw new IndexOutOfBoundsException( "Specified length causes crossing of origin: " + length + "; " + toString() );
		}
	}
	
	/**
	 * Return the adjacent location of specified length directly downstream of this location.
	 *
	 * @return The downstream location.
	 * @param length The length of the downstream location.
	 * @throws IndexOutOfBoundsException Specified length causes crossing of origin.
	 */
	public Location downstream( int length )
	{
		if( length < 0 )
		{
			throw new IllegalArgumentException( "Parameter must be >= 0; is=" + length );
		}
		
		if( Math.signum( mEnd + length) == Math.signum( mEnd ) || 0 == Math.signum( mEnd + length ) )
		{
			return new Location( mEnd, mEnd + length );
		}
		else
		{
			throw new IndexOutOfBoundsException( "Specified length causes crossing of origin: " + length + "; " + toString() );
		}
		
	}
	

	
	/**
	*   Return distance between this location and the other location.
	* 
	*	Distance is defined only if both locations are on same strand.
	 * 
	 * @param other The location to compare.
	 * @return The integer distance. Returns -1 if they overlap; 0 if directly adjacent.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public int distance( Location other )
	{
		if( isSameStrand( other ))
		{
			if( overlaps( other ))
			{
				return -1;
			}
			else
			{
				return ( mEnd <= other.mStart )? (other.mStart - mEnd) : (mStart - other.mEnd);
			}
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Return percent overlap of two locations.
	 * 
	 * @param other The location to compare.
	 * @return 100.0 * intersection(other).length() / this.length()
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public double percentOverlap( Location other )
	{
		if( length() > 0 && overlaps( other ))
		{
			return 100.0 * (((double) intersection( other ).length()) / (double) length());
		}
		else
		{
			return 0;
		}
	}
	
	/**	
	 * Check if this location and other location overlap.
	 *
	 * @param other The location to compare.
	 * @return True if they overlap.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean overlaps( Location other )
	{
		if( isSameStrand( other ))
		{
			return !( mStart >= other.mEnd || mEnd <= other.mStart );
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Check if this location contains the other.
	 *
	 * @param other The location to compare.
	 * @return True if other is entirely contained by this location. 
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean contains( Location other )
	{
		if( isSameStrand( other ))
		{
			return ( mStart <= other.mStart && mEnd >= other.mEnd );
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	} 
	
		
	/**	
	 * Check if this location starts after the other location starts.
	 * The locations may overlap.
	 *
	 * @param other The location to compare.
	 * @return True if this starts after other.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean startsAfter( Location other )
	{
		if( isSameStrand( other ))
		{
			return mStart > other.mStart;
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Check if this location starts before other location starts.
     * The locations may overlap.
     *
	 * @param other The location to compare.
	 * @return True if this starts before other.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean startsBefore( Location other )
	{
		if( isSameStrand( other ))
		{
			return mStart < other.mStart;
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Check if this location ends after other location ends.
     * The locations may overlap.
     *
	 * @param other The location to compare.
	 * @return True if location ends after other.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean endsAfter( Location other )
	{
		if( isSameStrand( other ) )
		{
			return mEnd > other.mEnd;
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Check if this location ends before other location ends.
     * The locations may overlap.
     *
	 * @param other The location to compare.
	 * @return True if this ends before other.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean endsBefore( Location other )
	{
		if( isSameStrand( other ) )
		{
			return mEnd < other.mEnd;
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Check if this location is entirely after the other location (no overlap).
	 *
	 * @param other The location to compare.
	 * @return True if this is after other.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean isAfter( Location other )
	{
		if( isSameStrand( other ) )
		{
			return mStart >= other.mEnd;
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}
	
	/**	
	 * Check if this location is entirely before other location (no overlap).
	 *
	 * @param other The location to compare.
	 * @return True if this is before other.
	 * @throws IllegalArgumentException Locations are on opposite strands.
	 */
	public boolean isBefore( Location other )
	{
		if( isSameStrand( other ) )
		{
			return mEnd <= other.mStart;
		}
		else
		{
			throw new IllegalArgumentException( "Locations are on opposite strands." );
		}
	}

	/**	
	 * Check if location is on negative strand.
	 * Note that Location( 0, 0 ) is by construction defined to be on the
	 * positive strand.
	 *
	 * @return True if on negative (reverse) strand.
	 */
	public boolean isNegative()
	{
		return ( mStart <= 0 && mEnd <= 0 );
	}
	
	/**	
	 * Return location that is in same position on opposite strand.
	 *
	 * @return Location on opposite strand.
	 */
	public Location opposite()
	{
		return new Location( - mEnd, - mStart );
	}
	
	/**	
	 * Check if this location is on same strand as other location.
	 *
	 * @param other The location to compare.
	 * @return True if on same strand.
	 */
	public boolean isSameStrand( Location other )
	{
		return ( isNegative() && other.isNegative() ) || ( !isNegative() && !other.isNegative() );
	}
	
		
	/**	
	 * Return a string representation of location.
	 *
	 * @return Text string.
	 */
	public String toString()
	{
		return new String( "[L=" + (mEnd - mStart) + "; S=" + mStart + "; E=" + mEnd +"]" );
	}
	
	/**	
	 * Compare locations for equality.
	 *
	 * @param other The location to compare.
	 * @return True if locations are the same.
	 */
	public boolean equals( Location other )
	{
		return mStart == other.mStart && mEnd == other.mEnd;
	}
	
	
	/**	
	 * 
	 */
	private boolean isHealthy()
	{
		return ( mStart <= mEnd ) && (( mStart <= 0 && mEnd <= 0 ) || (mStart >= 0 && mEnd >= 0));
	}
	
	//shorthand for testing
	static private Location L( int s, int e )
	{
	   return new Location( s, e );
	}
    
    static private Location R( int s, int e )
    {
       return new Location( -e, -s );
    }
    
    
	/**
	 * @deprecated
	 */
	public static void main( String args[] )
	throws Exception
	{
		Location p3_7= new Location( 3, 7 );
		Location p16_19= new Location( 16, 19 );
		Location p15_19= new Location( 15, 19 );
		Location p15_16= new Location( 15, 16 );
		Location p10_17= new Location( 10, 17 );
		Location p10_12= new Location( 10, 12 );
		Location p14_17= new Location( 14, 17 );
		Location p14_14= new Location( 14, 14 );
		
		Location r13_17= new Location( 13, 17 );
		Location r21_25= new Location( 21, 25 );
		
		Location r4_7= new Location( 4, 7 );
		Location r2_5= new Location( 2, 5 );
		Location r0_3= new Location( 0, 3 );
		Location r5_8= new Location( 5, 8 );
		
		//distance
		assert L(14,14).distance( L(3,7) ) == 7;
		assert L(3,7).distance( L(14,14) ) == 7;
		assert L(1,4).distance( L(7, 10) ) == 3;
		
		//union
		assert p10_12.union( p14_17 ).equals( p10_17 );
		assert p14_17.union( p10_12 ).equals( p10_17 );		
		assert p15_19.union( p15_16).equals( p15_19 );
		
		//intersection
		assert r13_17.union( r21_25 ).intersection( r21_25 ).equals( new Location( 21, 25 ));
		
				
		//isBefore
		assert r2_5.isBefore( r5_8 );		
		assert !r2_5.isBefore( r4_7 );
		
		//isAfter
		assert r5_8.isAfter( r2_5 );
		assert !r5_8.isAfter( r4_7 );
		
		//contains
		assert p15_19.contains( p16_19 );
				
		//overlaps		
		assert r2_5.overlaps( r4_7 );
		assert r2_5.overlaps( r0_3 );
		assert !r5_8.overlaps( r2_5 );
		assert !r2_5.overlaps( r5_8 );
		
		
        //prefix
        assert L(2,20).prefix(1).equals( L(2,3));
        assert L(2,20).prefix(-1).equals( L(2,19));
        assert L(2,20).prefix( L(10,12)).equals( L(2,10));
       
        //suffix
        assert L(2,20).suffix(1).equals( L(3,20));
        assert L(2,20).suffix(-1).equals( L(19,20));
        assert L(2,20).suffix( L(10,12)).equals( L(12,20));
       
        
        //upstream
        //downstream
        
        //startsBefore
        //startsAfter
        //endsBefore
        //endsAfter
        
        //equals
		
		//percentoverlap
		
		
		//plus
		//minus
		//isNegative
		//oppposite
		
		//fromBio, etc.
		
		System.out.println( "JavaGene.Location Passed." );
		
	}
	
}		
		 