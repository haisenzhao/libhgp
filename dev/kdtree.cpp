//kd_tree.h  kd_tree鐨勫ご鏂囦欢
//#include "StdAfx.h"

//澶存枃浠�
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#ifdef USE_LIST_NODE_ALLOCATOR

#ifndef NO_PTHREADS
#include <pthread.h>
#else

#ifndef I_WANT_THREAD_BUGS
#error "You are compiling with the fast list node allocator, with pthreads disabled! This WILL break if used from multiple threads."
#endif	/* I want thread bugs */

#endif	/* pthread support */
#endif	/* use list node allocator */


//瓒呭钩闈㈢殑缁撴瀯浣�
//鍖呮嫭涓€涓�灞炴€х殑缁存暟鍜屾瘡缁村潗鏍囩殑鏈€澶у拰鏈€灏忓€兼瀯鎴愮殑鏁扮粍
struct kdhyperrect {
	int dim;
	double *min, *max;              /* minimum/maximum coords */
};

//鑺傜偣鐨勭粨鏋勪綋锛屼篃灏辨槸浜嬩緥鐨勭粨鏋勪綋
struct kdnode {
	double *pos;
	int dir;
	void *data;

	struct kdnode *left, *right;	/* negative/positive side */
};

//杩斿洖缁撴灉鑺傜偣锛� 鍖呮嫭鏍戠殑鑺傜偣,璺濈�诲€�, 鏄�涓€涓�鍗曢摼琛ㄧ殑褰㈠紡
struct res_node {
	struct kdnode *item;
	double dist_sq;
	struct res_node *next;
};

//鏍戞湁鍑犱釜灞炴€э紝涓€鏄�缁存暟锛屼竴鏄�鏍戞牴鑺傜偣锛屼竴鏄�瓒呭钩闈�锛屼竴鏄�閿€姣乨ata鐨勫嚱鏁�
struct kdtree {
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
	void (*destr)(void*);
};

//kdtree鐨勮繑鍥炵粨鏋滐紝鍖呮嫭kdtree锛岃繖鏄�涓€涓�鍙岄摼琛ㄧ殑褰㈠紡
struct kdres {
	struct kdtree *tree;
	struct res_node *rlist, *riter;  //鍙岄摼琛�?
	int size;
};

//璁＄畻骞虫柟鐨勫畯瀹氫箟,鐩稿綋浜庡嚱鏁�
#define SQ(x)			((x) * (x))


static void clear_rec(struct kdnode *node, void (*destr)(void*));
static int insert_rec(struct kdnode **node, const double *pos, void *data, int dir, int dim);
static int rlist_insert(struct res_node *list, struct kdnode *item, double dist_sq);
static void clear_results(struct kdres *set);

static struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max);
static void hyperrect_free(struct kdhyperrect *rect);
static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect);
static void hyperrect_extend(struct kdhyperrect *rect, const double *pos);
static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos);

#ifdef USE_LIST_NODE_ALLOCATOR
static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node*);
#else
#define alloc_resnode()		malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif


//鍒涘缓涓€涓猭dtree
struct kdtree *kd_create(int k)
{
	struct kdtree *tree;

	if(!(tree = (kdtree*)malloc(sizeof *tree))) {
		return 0;
	}

	tree->dim = k;
	tree->root = 0;
	tree->destr = 0;
	tree->rect = 0;

	return tree;
}

//閲婃斁鎺塳dtree
void kd_free(struct kdtree *tree)
{
	if(tree) {
		kd_clear(tree);
		free(tree);
	}
}

//娓呴櫎鎺夎秴骞抽潰,鏄�鎸夎妭鐐归€掑綊鍦拌繘琛岀殑
static void clear_rec(struct kdnode *node, void (*destr)(void*))
{
	if(!node) return;   //涓€涓�鑺傜偣瀵瑰簲涓€涓�瓒呭钩闈�

	//閫掑綊鍑芥暟锛岄€掑綊鍦版竻闄ゆ帀浜屽弶鏍戝乏鍒嗘敮鐨勮秴骞抽潰鍜屼簩鍙夋爲鍙冲垎鏀�鐨勮秴骞抽潰
	clear_rec(node->left, destr);
	clear_rec(node->right, destr);
	
	//濡傛灉data娓呮�氬嚱鏁颁笉涓虹┖,灏遍噴鏀炬帀data
	if(destr) 
	{
		destr(node->data);
	}
	//閲婃斁鑺傜偣鐨勫潗鏍囨暟缁�
	free(node->pos);
	//閲婃斁鑺傜偣
	free(node);
}

//kdtree娓呴櫎
void kd_clear(struct kdtree *tree)
{
	//娓呴櫎鏍戜腑姣忎釜鑺傜偣鐨勮秴骞抽潰,閲婃斁鏍戜腑鐨勫悇涓�鑺傜偣
	clear_rec(tree->root, tree->destr);
	tree->root = 0;

	//濡傛灉鏍戠殑瓒呭钩闈㈡寚閽堜笉涓虹┖,瀵瑰叾杩涜�岄噴鏀�
	if (tree->rect) 
	{
		hyperrect_free(tree->rect);
		tree->rect = 0;
	}
}

//鏁版嵁閿€姣侊紝鐢ㄤ竴涓�澶栨潵鐨勫嚱鏁版潵杩涜�宒ata鐨勯攢姣�
void kd_data_destructor(struct kdtree *tree, void (*destr)(void*))
{
	//鐢ㄥ�栨潵鐨勫嚱鏁版潵鎵ц�宬dtree鐨勯攢姣佸嚱鏁�
	tree->destr = destr;
}


//鍦ㄤ竴涓�鏍戣妭鐐逛綅缃�澶勬彃鍏ヨ秴鐭╁舰
static int insert_rec(struct kdnode **nptr, const double *pos, void *data, int dir, int dim)
{
	int new_dir;
	struct kdnode *node;

	//濡傛灉杩欎釜鑺傜偣鏄�涓嶅瓨鍦ㄧ殑
	if(!*nptr) 
	{
		//鍒嗛厤涓€涓�缁撶偣
		if(!(node = (kdnode *)malloc(sizeof *node))) 
		{
			return -1;
		}
		if(!(node->pos = (double*)malloc(dim * sizeof *node->pos))) {
			free(node);
			return -1;
		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
		node->data = data;
		node->dir = dir;
		node->left = node->right = 0;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if(pos[node->dir] < node->pos[node->dir]) {
		return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
	}
	return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

//鑺傜偣鎻掑叆鎿嶄綔
//鍙傛暟涓�:瑕佽繘琛屾彃鍏ユ搷浣滅殑kdtree,瑕佹彃鍏ョ殑鑺傜偣鍧愭爣,瑕佹彃鍏ョ殑鑺傜偣鐨勬暟鎹�
int kd_insert(struct kdtree *tree, const double *pos, void *data)
{
	//鎻掑叆瓒呯煩褰�
	if (insert_rec(&tree->root, pos, data, 0, tree->dim)) 
	{
		return -1;
	}
	//濡傛灉鏍戣繕娌℃湁瓒呯煩褰�,灏卞垱寤轰竴涓�瓒呯煩褰�
	//濡傛灉宸茬粡鏈変簡瓒呯煩褰�,灏辨墿灞曞師鏈夌殑瓒呯煩褰�
	if (tree->rect == 0) 
	{
		tree->rect = hyperrect_create(tree->dim, pos, pos);
	} 
	else 
	{
		hyperrect_extend(tree->rect, pos);
	}

	return 0;
}

//鎻掑叆float鍨嬪潗鏍囩殑鑺傜偣
//鍙傛暟涓�:瑕佽繘琛屾彃鍏ユ搷浣滅殑kdtree,瑕佹彃鍏ョ殑鑺傜偣鍧愭爣,瑕佹彃鍏ョ殑鑺傜偣鐨勬暟鎹�
//灏唂loat鍨嬬殑鍧愭爣璧嬪€肩粰double鍨嬬殑缂撳啿鍖�,缁忚繃杩欎釜绫诲瀷杞�鍖栧悗杩涜�屾彃鍏�
//鏈�璐ㄤ笂鏄�涓€绉嶇被鍨嬭浆鍖�
int kd_insertf(struct kdtree *tree, const float *pos, void *data)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int res, dim = tree->dim;

	//濡傛灉kdtree鐨勭淮鏁板ぇ浜�16, 鍒嗛厤dim缁磀ouble绫诲瀷鐨勬暟缁�
	if(dim > 16) 
	{
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = (double*)alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (double*)malloc(dim * sizeof *bptr))) 
			{
				return -1;
			}
	} 
	//濡傛灉kdtree鐨勭淮鏁板皬浜�16, 鐩存帴灏嗘寚閽堟寚鍚戝凡鍒嗛厤鐨勫唴瀛�
	else 
	{
		bptr = buf = sbuf;
	}

	//灏嗚�佹彃鍏ョ偣鐨勪綅缃�鍧愭爣璧嬪€肩粰鍒嗛厤鐨勬暟缁�
	while(dim-- > 0) 
	{
		*bptr++ = *pos++;
	}

	//璋冪敤鑺傜偣鎻掑叆鍑芥暟kd_insert
	res = kd_insert(tree, buf, data);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
        //閲婃斁缂撳瓨
		free(buf);
	return res;
}

//缁欏嚭涓夌淮鍧愭爣鍊肩殑涓夌淮kdtree鎻掑叆
int kd_insert3(struct kdtree *tree, double x, double y, double z, void *data)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

//缁欏嚭涓夌淮float鍨嬪潗鏍囧€肩殑涓夌淮kdtree鎻掑叆
int kd_insert3f(struct kdtree *tree, float x, float y, float z, void *data)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_insert(tree, buf, data);
}

//鎵惧埌鏈€杩戦偦鐨勭偣
//鍙傛暟涓�:鏍戣妭鐐规寚閽�, 浣嶇疆鍧愭爣, 闃堝€�, 杩斿洖缁撴灉鐨勮妭鐐�, bool鍨嬫帓搴�,缁村害
static int find_nearest(struct kdnode *node, const double *pos, double range, struct res_node *list, int ordered, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;

	if(!node) return 0;  //娉ㄦ剰杩欎釜鍦版柟,褰撹妭鐐逛负绌虹殑鏃跺€�,琛ㄦ槑宸茬粡鏌ユ壘鍒版渶缁堢殑鍙跺瓙缁撶偣,杩斿洖鍊间负闆�

	dist_sq = 0;
	//璁＄畻涓や釜鑺傜偣闂寸殑骞虫柟鍜�
	for(i=0; i<dim; i++) 
	{
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	//濡傛灉璺濈�诲湪闃堝€艰寖鍥村唴,灏卞皢鍏舵彃鍏ュ埌杩斿洖缁撴灉閾捐〃涓�
	if(dist_sq <= SQ(range)) 
	{		
		if(rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1) 
		{
			return -1;
		}
		added_res = 1;
	}

	//鍦ㄨ繖涓�鑺傜偣鐨勫垝鍒嗘柟鍚戜笂,姹備袱鑰呬箣闂寸殑宸�鍊�
	dx = pos[node->dir] - node->pos[node->dir];

	//鏍规嵁杩欎釜宸�鍊肩殑绗﹀彿, 閫夋嫨杩涜�岄€掑綊鏌ユ壘鐨勫垎鏀�鏂瑰悜
	ret = find_nearest(dx <= 0.0 ? node->left : node->right, pos, range, list, ordered, dim);
	//濡傛灉杩斿洖鐨勫€煎ぇ浜庣瓑浜庨浂,琛ㄦ槑鍦ㄨ繖涓�鍒嗘敮涓�鏈夋弧瓒虫潯浠剁殑鑺傜偣,鍒欒繑鍥炵粨鏋滅殑涓�鏁拌繘琛岀疮鍔�,骞跺湪鑺傜偣鐨勫彟涓€涓�鏂瑰悜杩涜�屾煡鎵炬渶杩戠殑鑺傜偣
	if(ret >= 0 && fabs(dx) < range) 
	{
		added_res += ret;
		ret = find_nearest(dx <= 0.0 ? node->right : node->left, pos, range, list, ordered, dim);
	}
	if(ret == -1) 
	{
		return -1;
	}
	added_res += ret;

	return added_res;
}


//鎵惧埌鏈€杩戦偦鐨刵涓�鑺傜偣

#if 0
static int find_nearest_n(struct kdnode *node, const double *pos, double range, int num, struct rheap *heap, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;

	if(!node) return 0;
	
	/* if the photon is close enough, add it to the result heap */
	//濡傛灉瓒冲�熻繎灏卞皢鍏跺姞鍏ュ埌缁撴灉鍫嗕腑
	dist_sq = 0;
	//璁＄畻涓よ€呴棿鐨勬�у紡璺濈��
	for(i=0; i<dim; i++) 
	{
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	//濡傛灉璁＄畻鎵€寰楄窛绂诲皬浜庨槇鍊�
	if(dist_sq <= range_sq) {
	//濡傛灉鍫嗙殑澶у皬澶т簬num,涔熷氨鏄�澶т簬鎬荤殑瑕佹壘鐨勮妭鐐规暟
		if(heap->size >= num)
		{
			/* get furthest element */
			//寰楀埌鏈€杩滅殑鑺傜偣
			struct res_node *maxelem = rheap_get_max(heap);

			/* and check if the new one is closer than that */
			//娴嬭瘯杩欎釜鑺傜偣鏄�涓嶆槸姣旀渶杩滅殑鑺傜偣瑕佽繎
			if(maxelem->dist_sq > dist_sq) 
			{
			//濡傛灉鏄�鐨勮瘽,灏辩Щ闄ゆ渶杩滅殑鑺傜偣
				rheap_remove_max(heap);
				//骞跺皢姝よ妭鐐规彃鍏ュ爢涓�
				if(rheap_insert(heap, node, dist_sq) == -1) 
				{
					return -1;
				}
				added_res = 1;

				range_sq = dist_sq;
			}
		} 
		//濡傛灉鍫嗙殑澶у皬灏忎簬num,鐩存帴灏嗘�よ妭鐐规彃鍏ュ爢涓�
		else 
		{
			if(rheap_insert(heap, node, dist_sq) == -1) 
			{
				return =1;
			}
			added_res = 1;
		}
	}


	/* find signed distance from the splitting plane */
	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest_n(dx <= 0.0 ? node->left : node->right, pos, range, num, heap, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest_n(dx <= 0.0 ? node->right : node->left, pos, range, num, heap, dim);
	}
}
#endif


static void kd_nearest_i(struct kdnode *node, const double *pos, struct kdnode **result, double *result_dist_sq, struct kdhyperrect* rect)
{
	int dir = node->dir;
	int i;
	double dummy, dist_sq;
	struct kdnode *nearer_subtree, *farther_subtree;
	double *nearer_hyperrect_coord, *farther_hyperrect_coord;

	/* Decide whether to go left or right in the tree */
	//鍦ㄤ簩鍙夋爲涓�,鍐冲畾鍚戝乏璧拌繕鏄�鍚戝彸璧�
	dummy = pos[dir] - node->pos[dir];
	if (dummy <= 0) 
	{
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = rect->max + dir;
		farther_hyperrect_coord = rect->min + dir;
	} 
	else 
	{
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = rect->min + dir;
		farther_hyperrect_coord = rect->max + dir;
	}

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		kd_nearest_i(nearer_subtree, pos, result, result_dist_sq, rect);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

	/* Check the distance of the point at the current node, compare it
	 * with our best so far */
	dist_sq = 0;
	for(i=0; i < rect->dim; i++) 
	{
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if (dist_sq < *result_dist_sq) 
	{
		*result = node;
		*result_dist_sq = dist_sq;
	}

	if (farther_subtree) {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our
		 * minimum distance in result_dist_sq. */
		if (hyperrect_dist_sq(rect, pos) < *result_dist_sq) {
			/* Recurse down into farther subtree */
			kd_nearest_i(farther_subtree, pos, result, result_dist_sq, rect);
		}
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
	}
}

//姹俴dtree涓�涓庣偣pos鏈€杩戦偦鐨勫€�
struct kdres *kd_nearest(struct kdtree *kd, const double *pos)
{
	struct kdhyperrect *rect;
	struct kdnode *result;
	struct kdres *rset;
	double dist_sq;
	int i;

	//濡傛灉kd涓嶅瓨鍦�,鎴栬€呭叾瓒呭钩闈�涓嶅瓨鍦ㄧ殑璇�,鍒欏氨涓嶄細鏈夌粨鏋�
	if (!kd) return 0;
	if (!kd->rect) return 0;

	/* Allocate result set */
	//涓鸿繑鍥炵粨鏋滈泦鍚堝垎閰嶇┖闂�
	if(!(rset = (kdres*)malloc(sizeof *rset))) 
	{
		return 0;
	}
	if(!(rset->rlist = (res_node*)alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	/* Duplicate the bounding hyperrectangle, we will work on the copy */
	//澶嶅埗杈圭晫瓒呭钩闈�
	if (!(rect = hyperrect_duplicate(kd->rect))) 
	{
		kd_res_free(rset);
		return 0;
	}

	/* Our first guesstimate is the root node */
	result = kd->root;
	dist_sq = 0;
	for (i = 0; i < kd->dim; i++)
		dist_sq += SQ(result->pos[i] - pos[i]);

	/* Search for the nearest neighbour recursively */
	//閫掑綊鍦版煡鎵炬渶杩戦偦鐨勯偦灞�
	kd_nearest_i(kd->root, pos, &result, &dist_sq, rect);

	/* Free the copy of the hyperrect */
	//閲婃斁瓒呯煩褰�
	hyperrect_free(rect);

	/* Store the result */
	//瀛樺偍缁撴灉
	if (result) 
	{
		if (rlist_insert(rset->rlist, result, -1.0) == -1) 
		{
			kd_res_free(rset);
			return 0;
		}
		rset->size = 1;
		kd_res_rewind(rset);
		return rset;
	} 
	else 
	{
		kd_res_free(rset);
		return 0;
	}
}

//kd_nearest鐨刦loat鐗逛緥
struct kdres *kd_nearestf(struct kdtree *tree, const float *pos)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int dim = tree->dim;
	struct kdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = (double*)alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (double*)malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_nearest(tree, buf);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
		free(buf);
	return res;
}

//kd_nearest鐨勪笁鍧愭爣鐗逛緥
struct kdres *kd_nearest3(struct kdtree *tree, double x, double y, double z)
{
	double pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return kd_nearest(tree, pos);
}

//kd_nearest鐨勪笁鍧愭爣float鐗逛緥
struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z)
{
	double pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return kd_nearest(tree, pos);
}

/* ---- nearest N search ---- */
/*
static kdres *kd_nearest_n(struct kdtree *kd, const double *pos, int num)
{
	int ret;
	struct kdres *rset;

	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	if((ret = find_nearest_n(kd->root, pos, range, num, rset->rlist, kd->dim)) == -1) {
		kd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	kd_res_rewind(rset);
	return rset;
}*/

//鎵惧埌婊¤冻璺濈�诲皬浜巖ange鍊肩殑鑺傜偣
struct kdres *kd_nearest_range(struct kdtree *kd, const double *pos, double range)
{
	int ret;
	struct kdres *rset;

	if(!(rset = (kdres*)malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = (res_node*)alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	if((ret = find_nearest(kd->root, pos, range, rset->rlist, 0, kd->dim)) == -1) {
		kd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	kd_res_rewind(rset);
	return rset;
}

//kd_nearest_range鐨刦loat鐗逛緥
struct kdres *kd_nearest_rangef(struct kdtree *kd, const float *pos, float range)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int dim = kd->dim;
	struct kdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = (double*)alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = (double*)malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = kd_nearest_range(kd, buf, range);
#ifndef NO_ALLOCA
	if(kd->dim > 256)
#else
	if(kd->dim > 16)
#endif
		free(buf);
	return res;
}

//kd_nearest_range鐨勪笁鍧愭爣鐗逛緥
struct kdres *kd_nearest_range3(struct kdtree *tree, double x, double y, double z, double range)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}

//kd_nearest_range鐨勪笁鍧愭爣float鐗逛緥
struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return kd_nearest_range(tree, buf, range);
}

//杩斿洖缁撴灉鐨勯噴鏀�
void kd_res_free(struct kdres *rset)
{
	clear_results(rset);
	free_resnode(rset->rlist);
	free(rset);
}

//鑾峰彇杩斿洖缁撴灉闆嗗悎鐨勫ぇ灏�
int kd_res_size(struct kdres *set)
{
	return (set->size);
}

//鍐嶆�″洖鍒拌繖涓�鑺傜偣鏈�韬�鐨勪綅缃�
void kd_res_rewind(struct kdres *rset)
{
	rset->riter = rset->rlist->next;
}

//鎵惧埌杩斿洖缁撴灉涓�鐨勬渶缁堣妭鐐�
int kd_res_end(struct kdres *rset)
{
	return rset->riter == 0;
}

//杩斿洖缁撴灉鍒楄〃涓�鐨勪笅涓€涓�鑺傜偣
int kd_res_next(struct kdres *rset)
{
	rset->riter = rset->riter->next;
	return rset->riter != 0;
}

//灏嗚繑鍥炵粨鏋滅殑鑺傜偣鐨勫潗鏍囧拰data鎶藉彇鍑烘潵
void *kd_res_item(struct kdres *rset, double *pos)
{
	if(rset->riter) {
		if(pos) {
			memcpy(pos, rset->riter->item->pos, rset->tree->dim * sizeof *pos);
		}
		return rset->riter->item->data;
	}
	return 0;
}

//灏嗚繑鍥炵粨鏋滅殑鑺傜偣鐨勫潗鏍囧拰data鎶藉彇鍑烘潵,鍧愭爣涓篺loat鍨嬬殑鍊�
void *kd_res_itemf(struct kdres *rset, float *pos)
{
	if(rset->riter) {
		if(pos) {
			int i;
			for(i=0; i<rset->tree->dim; i++) {
				pos[i] = rset->riter->item->pos[i];
			}
		}
		return rset->riter->item->data;
	}
	return 0;
}

//灏嗚繑鍥炵粨鏋滅殑鑺傜偣鐨勫潗鏍囧拰data鎶藉彇鍑烘潵,鍧愭爣鍏蜂綋褰㈠紡缁欏嚭
void *kd_res_item3(struct kdres *rset, double *x, double *y, double *z)
{
	if(rset->riter) {
		if(*x) *x = rset->riter->item->pos[0];
		if(*y) *y = rset->riter->item->pos[1];
		if(*z) *z = rset->riter->item->pos[2];
	}
	return 0;
}

//灏嗚繑鍥炵粨鏋滅殑鑺傜偣鐨勫潗鏍囧拰data鎶藉彇鍑烘潵,鍧愭爣涓篺loat鍨嬬殑鍊�,鍧愭爣鍏蜂綋褰㈠紡缁欏嚭
void *kd_res_item3f(struct kdres *rset, float *x, float *y, float *z)
{
	if(rset->riter) {
		if(*x) *x = rset->riter->item->pos[0];
		if(*y) *y = rset->riter->item->pos[1];
		if(*z) *z = rset->riter->item->pos[2];
	}
	return 0;
}

//鑾峰彇data鏁版嵁
void *kd_res_item_data(struct kdres *set)
{
	return kd_res_item(set, 0);
}

/* ---- hyperrectangle helpers ---- */
//鍒涘缓瓒呭钩闈�,鍖呮嫭涓変釜鍙傛暟:缁村害,姣忕淮鐨勬渶灏忓€煎拰鏈€澶у€兼暟缁�
static struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max)
{
	size_t size = dim * sizeof(double);
	struct kdhyperrect* rect = 0;

	if (!(rect = (kdhyperrect*)malloc(sizeof(struct kdhyperrect)))) 
	{
		return 0;
	}

	rect->dim = dim;
	if (!(rect->min = (double*)malloc(size))) {
		free(rect);
		return 0;
	}
	if (!(rect->max = (double*)malloc(size))) {
		free(rect->min);
		free(rect);
		return 0;
	}
	memcpy(rect->min, min, size);
	memcpy(rect->max, max, size);

	return rect;
}

//閲婃斁瓒呭钩闈㈢粨鏋勪綋
static void hyperrect_free(struct kdhyperrect *rect)
{
	free(rect->min);
	free(rect->max);
	free(rect);
}

//璧嬪€艰秴骞抽潰缁撴瀯浣�
static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect)
{
	return hyperrect_create(rect->dim, rect->min, rect->max);
}

//鏇存柊瓒呭钩闈㈢粨鏋勪綋鏈€澶�\鏈€灏忓€兼暟缁�
static void hyperrect_extend(struct kdhyperrect *rect, const double *pos)
{
	int i;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			rect->min[i] = pos[i];
		}
		if (pos[i] > rect->max[i]) {
			rect->max[i] = pos[i];
		}
	}
}

//璁＄畻鍥哄畾鍧愭爣鐐逛笌瓒呭钩闈�涔嬮棿鐨勮窛绂�
static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos)
{
	int i;
	double result = 0;

	for (i=0; i < rect->dim; i++) 
	{
		if (pos[i] < rect->min[i]) 
		{
			result += SQ(rect->min[i] - pos[i]);
		} 
		else if (pos[i] > rect->max[i]) 
		{
			result += SQ(rect->max[i] - pos[i]);
		}
	}
	return result;
}


/* ---- static helpers ---- */
#ifdef USE_LIST_NODE_ALLOCATOR
/* special list node allocators. */
static struct res_node *free_nodes;

#ifndef NO_PTHREADS
static pthread_mutex_t alloc_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

//鍒涘缓杩斿洖缁撴灉鑺傜偣
static struct res_node *alloc_resnode(void)
{
	struct res_node *node;

#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	if(!free_nodes) {
		node = malloc(sizeof *node);
	} else {
		node = free_nodes;
		free_nodes = free_nodes->next;
		node->next = 0;
	}

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif

	return node;
}

//閲婃斁杩斿洖缁撴灉鑺傜偣
static void free_resnode(struct res_node *node)
{
#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	node->next = free_nodes;
	free_nodes = node;

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif
}
#endif	/* list node allocator or not */


/* inserts the item. if dist_sq is >= 0, then do an ordered insert */
/* TODO make the ordering code use heapsort */
//鍑芥暟鍙傛暟: 杩斿洖缁撴灉鑺傜偣鎸囬拡,鏍戣妭鐐规寚閽�,璺濈�诲嚱鏁�
//灏嗕竴涓�缁撴灉鑺傜偣鎻掑叆鍒拌繑鍥炵粨鏋滅殑鍒楄〃涓�
static int rlist_insert(struct res_node *list, struct kdnode *item, double dist_sq)
{
	struct res_node *rnode;

	//鍒涘缓涓€涓�杩斿洖缁撴灉鐨勮妭鐐�
	if(!(rnode = (res_node*)alloc_resnode())) 
	{
		return -1;
	}
	rnode->item = item;           //瀵瑰簲鐨勬爲鑺傜偣
	rnode->dist_sq = dist_sq;     //瀵瑰簲鐨勮窛绂诲€�

	//褰撹窛绂诲ぇ浜庨浂鐨勬椂鍊�
	if(dist_sq >= 0.0) 
	{
		while(list->next && list->next->dist_sq < dist_sq) 
		{
			list = list->next;
		}
	}
	rnode->next = list->next;
	list->next = rnode;
	return 0;
}

//娓呴櫎杩斿洖缁撴灉鐨勯泦鍚�
//鏈�璐ㄤ笂鏄�涓�鍙岄摼琛ㄤ腑鍗曢摼琛ㄧ殑娓呯悊
static void clear_results(struct kdres *rset)
{
	struct res_node *tmp, *node = rset->rlist->next;

	while(node) 
	{
		tmp = node;
		node = node->next;
		free_resnode(tmp);
	}

	rset->rlist->next = 0;
}
