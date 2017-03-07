/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.genome.parsers.gff;

import org.biojava.nbio.genome.App;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Move a sliding window over a Location.
 * Window size and increment can be specified.
 * If the increment is negative, the iteration starts
 * at end of Location and moves toward beginning.
 *
 * @author Hanno Hinsch
 */
public class LocIterator implements Iterator<Location> {
	private static final Logger logger = LoggerFactory.getLogger(App.class);

	Location mBounds;
	int mPosition;
	int mWindowSize;
	int mIncrement;

	@SuppressWarnings("unused")
	private LocIterator() {};

	/**
	 * Construct an iterator that slides a window over a Location.
	 *
	 * @param bounds The location over which to iterate.
	 * @param windowSize The size of the moving window.
	 * @param increment The increment by which to move the window at each iteration.
	 * If increment is positive, the iteration starts at the beginning of the bounding location
	 * and moves toward the end; if the increment is negative, the iteration starts at the end and moves
	 * toward the begnning.
	 */
	public LocIterator( Location bounds, int windowSize, int increment )
	{
		mWindowSize= windowSize;
		mIncrement= increment;
		mBounds= bounds;

		if( windowSize <= 0 )
		{
			throw new IllegalArgumentException( "Window size must be positive." );
		}

		if( increment == 0 )
		{
			throw new IllegalArgumentException( "Increment must be non-zero." );
		}

		mPosition= 0;

	}



	/**
	 * Check if next window of specified size is available.
	 *
	 * @param windowSize Size of window. May be smaller or larger than default window size.
	 * @param increment The increment by which to move the window at each iteration. Note that this
	 * method does not actually change the position. However, it checks the sign of the increment parameter to determine
	 * the direction of the iteration.
	 * @return True if window of specified size is available.
	 * @throws IllegalArgumentException Window size parameter was not positive.
	 */
	public boolean hasNext( int windowSize, int increment )
	{
		if( windowSize <= 0 )
		{
			throw new IllegalArgumentException( "Window size must be positive." );
		}

		try
		{
			if( increment > 0 )
			{
				return windowSize == mBounds.suffix( mPosition ).prefix( windowSize ).length();
			}
			else
			{
				if( mPosition == 0 )
				{
					return windowSize == mBounds.suffix( - windowSize ).length();
				}
				else
				{
					return windowSize == mBounds.prefix( mPosition ).suffix( - windowSize ).length();
				}
			}
		}
		catch( Exception e )
		{
			return false;
		}
	}

	/**
	 * Check if next window of default size is available.
	 *
	 * @return True if window of default size is available. The default size
	 * is the size specified in the LocIterator constructor.
	 */
	@Override
	public boolean hasNext()
	{
		return hasNext( mWindowSize, mIncrement );
	}

	/**
	 * Get portion of bounding location that has not yet been retrieved by next() method.
	 *
	 * @return The location not yet retrieved.
	 */
	public Location remainder()
	{
		Location remainder= null;

		if( mPosition == 0 )
		{
			remainder= mBounds;
		}
		else
		{
			if( mIncrement > 0 )
			{
				remainder = mBounds.suffix( mPosition );
			}
			else
			{
				remainder = mBounds.prefix( mPosition );
			}
		}

		return remainder;

	}

	/**
	 * Get next window of default size, then increment position by default amount. Both
	 * defaults are specified in the LocIterator constructor.
	 *
	 * @return Location of next window.
	 * @throws IndexOutOfBoundsException The next window was not within the bounding location.
	 */
	@Override
	public Location next()
	{
        if(!hasNext()){
            throw new NoSuchElementException();
        }
        return next( mWindowSize, mIncrement );
	}

	/**
	 * Get next window of specified size, then increment position by specified amount.
	 *
	 * @return Location of next window.
	 * @param windowSize Size of window to get.
	 * @param increment Amount by which to shift position. If increment is positive, the position is shifted
	 * toward the end of the bounding location; if increment is negative, the position is shifted toward
	 * the beginning of the bounding location.
	 * @throws IndexOutOfBoundsException The next window was not within the bounding location.
	 * @throws IllegalArgumentException The increment was zero, or windowSize was not positive.
	 */
	public Location next( int windowSize, int increment )
	{
		if( windowSize <= 0 )
		{
			throw new IllegalArgumentException( "Window size must be positive." );
		}

		if( increment == 0 )
		{
			throw new IllegalArgumentException( "Increment must be non-zero." );
		}

		Location r;

		try
		{
			if( increment > 0 )
			{
				r= mBounds.suffix( mPosition ).prefix( windowSize );
			}
			else
			{
				if( mPosition == 0 )
				{
					r= mBounds.suffix( - windowSize );
				}
				else
				{
					r= mBounds.prefix( mPosition ).suffix( - windowSize );
				}
			}

			mPosition+= increment;

		}
		catch( Exception e )
		{
			throw new IndexOutOfBoundsException( e.toString() );
		}

		return r;
	}

	/**
	 * Get string representation of iterator.
	 *
	 * @return Description of internal state.
	 */
	@Override
	public String toString()
	{
		return "bounds=" + mBounds.toString() + "; pos=" + mPosition + "; winsize=" + mWindowSize + "; inc=" + mIncrement;
	}

	/**
	 * Unsupported.
	 *
	 * @throws UnsupportedOperationException
	 */
	@Override
	public void remove()
	{
		throw new UnsupportedOperationException();
	}

	/**
	 * @deprecated
	 */
	@Deprecated
	public static void main(String[] args )
	{

		Location r= new Location( 10, 21 );

		logger.info( "10 to 21, 1 by 1" );
		for( Location t: r.window( 1, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 3" );
		for( Location t: r.window( 3, 3 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 1" );
		for( Location t: r.window( 3, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 11 by 1" );
		for( Location t: r.window( 11, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 12 by 1" );
		for( Location t: r.window( 12, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by -1" );
		for( Location t: r.window( 1, -1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by -3" );
		for( Location t: r.window( 3, -3 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by -1" );
		for( Location t: r.window( 3, -1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by 1" );
		for( Location t: r.window( 1, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 3" );
		for( Location t: r.window( 3, 3 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 1" );
		for( Location t: r.window( 3, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 11 by 1" );
		for( Location t: r.window( 11, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 12 by 1" );
		for( Location t: r.window( 12, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by -1" );
		for( Location t: r.window( 1, -1 )) { logger.info( t.toString() ); }

		//reverse strand
		r= r.opposite();
		logger.info( "reverse strand" );

		logger.info( "10 to 21, 1 by 1" );
		for( Location t: r.window( 1, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 3" );
		for( Location t: r.window( 3, 3 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 1" );
		for( Location t: r.window( 3, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 11 by 1" );
		for( Location t: r.window( 11, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 12 by 1" );
		for( Location t: r.window( 12, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by -1" );
		for( Location t: r.window( 1, -1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by -3" );
		for( Location t: r.window( 3, -3 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by -1" );
		for( Location t: r.window( 3, -1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by 1" );
		for( Location t: r.window( 1, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 3" );
		for( Location t: r.window( 3, 3 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 3 by 1" );
		for( Location t: r.window( 3, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 11 by 1" );
		for( Location t: r.window( 11, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 12 by 1" );
		for( Location t: r.window( 12, 1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by -1" );
		for( Location t: r.window( 1, -1 )) { logger.info( t.toString() ); }

		logger.info( "10 to 21, 1 by 1 (+2)" );
		LocIterator i= r.iterator();
		int chunk= 1;
		while( i.hasNext( 1, chunk ) )
		{
			Location t= i.next( 1, chunk );
			logger.info( t.toString() );
			chunk+= 2;
		}

		//FIXME test remainder()

		logger.info("JavaGene.LocIterator Passed.");
	}

}
