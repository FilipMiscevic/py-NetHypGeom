from rpy2.robjects.packages import importr
from random_walk import *
from itertools import permutations
from plotly.offline import plot
from plotly.graph_objs import *

base = importr('base')
stats = importr('stats')

### R wrapper functions ###
string = '''labne_hm_from_matrix <- function(net, gma = NA, Temp = 0.1, k.speedup = 10, m.in = NA, L.in = NA, w = "auto"){

  m=as.matrix(net>0) # coerces the data set as a matrix
  diag(m) <- 0
  g=graph.adjacency(m,mode='undirected',weighted=NULL)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                edge.attr.comb = igraph_opt('edge.attr.comb'))

  return(labne_hm(g, gma = gma, Temp = Temp, k.speedup = k.speedup, m.in = m.in, L.in = L.in, w = w))
  }
'''

pp = SignatureTranslatedAnonymousPackage(string,'pp')
### R wrapper functions ###


def plot_hyperbolic_network(mat,polar_coords,colors=None,labels=[]):
    '''
    Plot a network
    :param mat:
    :param polar_coords:
    :param colors:
    :param labels:
    :return:
    '''

    ut = np.triu(mat)
    edges = np.array( np.where(ut != 0)).T

    xs = polar_coords[:,0]*np.cos(polar_coords[:,1])
    ys = polar_coords[:,0]*np.sin(polar_coords[:,1])

    xys = zip(xs,ys)

    edge_trace = Scatter(
        x=[],
        y=[],
        line=Line(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    for e1,e2 in edges:
        x0, y0 = xys[e1]
        x1, y1 = xys[e2]
        edge_trace['x'] += [x0, x1, None]
        edge_trace['y'] += [y0, y1, None]

    node_trace = Scatter(
        x=[],
        y=[],
        text=labels,
        mode='markers',
        hoverinfo='text',
        marker=Marker(
            showscale=True,
            # colorscale options
            # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
            # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
            colorscale='YIGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2)))

    for (x,y) in xys:
        node_trace['x'].append(x)
        node_trace['y'].append(y)

    fig = Figure(data=Data([edge_trace, node_trace]),
                 layout=Layout(
                     title='<br>Network in the hyperbolic disc',
                     titlefont=dict(size=16),
                     showlegend=False,
                     hovermode='closest',
                     margin=dict(b=20, l=5, r=5, t=40),
                     annotations=[dict(
                         showarrow=False,
                         xref="paper", yref="paper",
                         x=0.005, y=-0.002)],
                     xaxis=XAxis(showgrid=False, zeroline=False, showticklabels=False),
                     yaxis=YAxis(showgrid=False, zeroline=False, showticklabels=False)))

    plot(fig, filename='networkx')



def greedy_route_packets(mat,dist,pairs=None):

    N = len(dist)

    if pairs is None:
        all_pairs = np.array(list(permutations(xrange(N),2)))
    else:
        all_pairs = pairs

    #all_pairs = all_pairs[np.random.choice(xrange(len(all_pairs)),10000,replace=False)]

    hops = np.zeros(len(all_pairs))

    j = 0
    for (source,target) in all_pairs:
        previous = -1

        while True:

            # assume there is strictly one closer neighbor
            neighbors = np.where( mat[source] != 0 )[0]
            closest = neighbors[np.argmin(dist[target,neighbors])]
            #print closest,source,target

            if closest == target:
                hops[j] += 1
                break
            elif closest == previous:
                hops[j] = 0
                break
            else:
                previous = source
                source = closest
                hops[j] += 1
        j += 1

    return hops


#########################
### DISTANCE MEASURES ###
#########################


def euclidean_distance(x1, x2=None):
    '''
    Compute distances in Euclidean space.
    :param x1: Rows of points in Euclidean space of arbitrary dimension.
    :param x2: Rows of points in Euclidean space to compute distance to. If x2 is None, x1 is used, thereby computing pairwise distances.
    :return: Distance matrix.
    '''
    if x2 is None:
        x2 = x1
    # numpy broadcasting trick to compute pairwise distances
    return np.sqrt((((x1[np.newaxis, :, :] - x2[:, np.newaxis, :]) ** 2).sum(2)).T)


def hyperbolic_distance(x1, x2=None):
    '''
    Compute distances in hyperbolic space.
    :param x1: Rows of points in the H^2 hyperbolic model. The first coordinate is radius and the second is theta, the angular coordinate.
    :param x2: Rows of points in the same format as x1 to compute distance to. If x2 is None, x1 is used, thereby computing pairwise distances.
    :return: Distance matrix.
    '''
    if x2 is None:
        x2 = x1
    # numpy broadcasting trick to compute pairwise distances
    #d = np.arccosh(np.cosh(x1[0])*np.cosh(x2[0]) - np.sinh(x1[0])*np.sinh(x2[0])*np.cos(np.pi-np.abs(np.pi-np.abs(x1[1]-x2[1]))))
    d = ((np.cosh(x1[:,0,np.newaxis][np.newaxis, :, :])*np.cosh(x2[:,0,np.newaxis][:, np.newaxis, :]) - np.sinh(x1[:,0,np.newaxis][np.newaxis, :, :])*np.sinh(x2[:,0,np.newaxis][:, np.newaxis, :])*np.cos(np.pi-np.abs(np.pi-np.abs(x1[:,1,np.newaxis][np.newaxis, :, :]-x2[:,1,np.newaxis][:, np.newaxis, :])))).sum(2)).T

    # Due to precision problems, distances less than 1 should be 1 to get correct hyperbolic distance of 0
    d[d<1] = 1

    d = np.arccosh(d)

    # distances to self should be 0
    np.fill_diagonal(d,0)

    return d


if __name__ == "__main__":
    # select subject
    s = 3

    d = DatabaseConnection(database='lau2')

    table = 'right_hemisphere_234'
    print table



    iq = []
    e_diff_norm = []
    e_rout_norm = []
    e_rout = []
    e_diff = []
    cc = []
    cc_norm = []
    pl = []
    pl_norm = []

    W_Ps = []

    ALL_DATA = {}

    #for table in tables:


    ids = []
    EDR = []
    HMR = []
    FLR = []
    scores = []
    hc = []

    q = 'SELECT `Subject ID`,`Original` from {} '.format(table, s)

    c, r = d.execute_and_return_query(q)
    if True:
        for i in range(len(r)):
            try:
                id = r[i][0]

                print id

                W_P = read_msgpack(str(r[i][1]))

                W_Ps.append( W_P )
                ids.append(id)



                ED = W_P.ED.as_matrix()
                W = W_P.W.as_matrix()
                np.fill_diagonal(W,0)
                np.fill_diagonal(ED,0)

                angle = pi/12


                hops = greedy_route_packets(W, ED)

                EDR.append(hops)

                print "ED routing: ", np.count_nonzero(hops) * 1.0 / len(hops)

                #test = pp.labne_hm_from_matrix(net=W, gma=2.3, Temp=0.75, k_speedup=0, w=angle)
                test = pp.labne_hm_from_matrix(net=W, gma=2.7, Temp=0.55, k_speedup=0, w=angle)
                test_df = pandas2ri.ri2py(test[1])
                hyperbolic_coords = test_df.as_matrix(['r', 'theta'])
                hc.append(hyperbolic_coords)

                gg = hyperbolic_distance(hyperbolic_coords, hyperbolic_coords)

                hops2 = greedy_route_packets(W, gg)

                HMR.append(hops2)

                print "HM routing: ", np.count_nonzero(hops2) * 1.0 / len(hops2)
                #plot_hyperbolic_network(W, test_df.as_matrix(['r','theta']), labels=W_P.major_axis)

                L = W_P.L.as_matrix()
                np.fill_diagonal(L,0)

                FL = interpolate_fibre_length(L,ED)

                hops3 = greedy_route_packets(W, FL)
                FLR.append(hops3)
                print "FL routing: ", np.count_nonzero(hops3) * 1.0 / len(hops3)

                scores.append([np.count_nonzero(hops) * 1.0 / len(hops),np.count_nonzero(hops3) * 1.0 / len(hops3),np.count_nonzero(hops2) * 1.0 / len(hops2)])

            except Exception as e:
                print (e, i, id)
                #print W_P
                continue

        print 'done!'
