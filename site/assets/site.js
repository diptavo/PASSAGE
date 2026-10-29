document.querySelectorAll('pre').forEach((block) => {
  if (!navigator.clipboard) return;
  const button = document.createElement('button');
  button.className = 'copy';
  button.type = 'button';
  button.textContent = 'Copy';
  button.setAttribute('aria-label', 'Copy code to clipboard');
  button.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(block.querySelector('code').textContent);
      button.textContent = 'Copied';
    } catch (_) {
      button.textContent = 'Select code to copy';
    }
    setTimeout(() => { button.textContent = 'Copy'; }, 2000);
  });
  block.appendChild(button);
});
