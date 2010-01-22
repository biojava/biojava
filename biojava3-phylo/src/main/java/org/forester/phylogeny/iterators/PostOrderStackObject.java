// $Id: PostOrderStackObject.java,v 1.6 2009/01/13 19:49:30 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny.iterators;

import org.forester.phylogeny.PhylogenyNode;

/*
 * @author Christian M. Zmasek
 * 
 * @version 1.00 -- last modified: 06/15/00
 */
public class PostOrderStackObject {

    final private PhylogenyNode _node;
    final private int           _phase;

    /**
     * Creates a new PostOrderStackObject object.
     * 
     * @param n
     *            DOCUMENT ME!
     * @param i
     *            DOCUMENT ME!
     */
    public PostOrderStackObject( final PhylogenyNode n, final int i ) {
        _node = n;
        _phase = i;
    }

    /**
     * DOCUMENT ME!
     * 
     * @return DOCUMENT ME!
     */
    public PhylogenyNode getNode() {
        return _node;
    }

    /**
     * DOCUMENT ME!
     * 
     * @return DOCUMENT ME!
     */
    public int getPhase() {
        return _phase;
    }
}
