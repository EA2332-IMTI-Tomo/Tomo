// a mutex object
pthread_mutex_t mutexCUDAInit;

//----------------------------------------------------------------------------------------------------------
void *run_func(void *arg)
{
  // calling cudaFree(0) within a critical section protected by the mutex
  pthread_mutex_lock(&mutexCUDAInit);
  cudaFree(0);
  pthread_mutex_unlock(&mutexCUDAInit);

  // your CUDA code here

  return NULL;
}

#define NUM_THREAD	10
//----------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  pthread_t pt[NUM_THREAD];
  int i;
  int status;

  // initializing a mutex
  pthread_mutex_init(&mutexCUDAInit, NULL);

  // create threads
  for(i = 0 ; i < NUM_THREAD ; i++){
    printf("create thread %d\n", i);
    status = pthread_create(&pt[i], NULL, run_func, NULL);
    if (status != 0) {
      perror("pthread_create");
      return 2;
    }
  }

  for(i = 0 ; i < NUM_THREAD ; i++){
    pthread_join(pt[i], NULL);
  }

  return 0;
}


